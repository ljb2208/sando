/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

// Lightweight modeling layer used by SolverHighsShot. It mirrors the subset of the Gurobi C++ API
// that SolverGurobi relies on (GRBVar, GRBLinExpr, GRBQuadExpr, GRBModel::addVar/addConstr/
// addQConstr/addGenConstrIndicator/remove/setObjective/optimize), so the trajectory formulation
// can be expressed identically. At optimize() time the model is dispatched to:
//   - HiGHS directly, for continuous LP/QP and for MILP, and
//   - SHOT (with HiGHS as its MIP subsolver), for convex MIQP and problems with quadratic
//     constraints.
// Indicator constraints, which neither HiGHS nor SHOT support natively, are reformulated with
// big-M constants computed from feasibility-based bound tightening on the linear constraints.
// After SHOT solves a mixed-integer problem, the binaries are fixed and the continuous problem is
// solved again (HiGHS for QPs), so the returned solution is not affected by big-M round-off.
// No third-party headers are included here; HiGHS and SHOT are only used in the .cpp file.

#pragma once

#include <atomic>
#include <chrono>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace hs {

constexpr double kInf = std::numeric_limits<double>::infinity();

// Variable types (same character codes as GRB_CONTINUOUS / GRB_BINARY)
constexpr char kContinuous = 'C';
constexpr char kBinary = 'B';

// Constraint senses (same character codes as GRB_LESS_EQUAL / GRB_GREATER_EQUAL / GRB_EQUAL)
constexpr char kLessEqual = '<';
constexpr char kGreaterEqual = '>';
constexpr char kEqual = '=';

// Objective senses (same values as GRB_MINIMIZE / GRB_MAXIMIZE)
constexpr int kMinimize = 1;
constexpr int kMaximize = -1;

/** @brief Optimization status, mirroring the Gurobi status codes used by SolverGurobi. */
enum class Status {
  kLoaded,       // optimize() has not been called since the last solve
  kOptimal,      // optimal (within the MIP gap for mixed-integer problems)
  kInfeasible,
  kInfOrUnbd,
  kUnbounded,
  kTimeLimit,
  kInterrupted,  // stopped by the termination flag
  kNumeric,
  kError,
};

/** @brief Which backend solved the last model. */
enum class Backend { kNone, kHighs, kShot };

/** @brief Exception thrown on modeling or solver errors (the counterpart of GRBException). */
class SolverException : public std::runtime_error {
 public:
  explicit SolverException(const std::string& msg) : std::runtime_error(msg) {}
};

/** @brief Handle to a decision variable of a Model. */
class Var {
 public:
  Var() = default;
  explicit Var(int index) : index_(index) {}
  int index() const { return index_; }
  bool valid() const { return index_ >= 0; }

 private:
  int index_{-1};
};

/** @brief Affine expression sum_i coeff_i * var_i + constant (the counterpart of GRBLinExpr). */
class LinExpr {
 public:
  LinExpr(double constant = 0.0) : constant_(constant) {}  // NOLINT: implicit like GRBLinExpr
  LinExpr(const Var& var, double coeff = 1.0);             // NOLINT: implicit like GRBLinExpr

  LinExpr& operator+=(const LinExpr& other);
  LinExpr& operator-=(const LinExpr& other);
  LinExpr& operator*=(double factor);
  LinExpr& operator/=(double divisor);

  /** @brief Number of (possibly repeated) variable terms. */
  size_t size() const { return vars_.size(); }
  int getVarIndex(size_t i) const { return vars_[i]; }
  double getCoeff(size_t i) const { return coeffs_[i]; }
  double getConstant() const { return constant_; }

 private:
  std::vector<int> vars_;
  std::vector<double> coeffs_;
  double constant_{0.0};
};

LinExpr operator+(LinExpr lhs, const LinExpr& rhs);
LinExpr operator-(LinExpr lhs, const LinExpr& rhs);
// Scalar overloads, so mixed expressions do not depend on implicit conversions
LinExpr operator+(LinExpr lhs, double rhs);
LinExpr operator+(double lhs, LinExpr rhs);
LinExpr operator-(LinExpr lhs, double rhs);
LinExpr operator-(double lhs, LinExpr rhs);
LinExpr operator-(LinExpr expr);
LinExpr operator*(LinExpr expr, double factor);
LinExpr operator*(double factor, LinExpr expr);
LinExpr operator/(LinExpr expr, double divisor);

/**
 * @brief Quadratic expression (cf. GRBQuadExpr):
 *   affine part + sum_k q_k * x_i(k) * x_j(k) + sum_m c_m * L_m * R_m,
 * where L_m, R_m are affine expressions. Products of affine expressions are kept unexpanded so
 * that sums of squares (e.g. ||jerk||^2) can be passed to SHOT in a well-scaled separable form.
 */
class QuadExpr {
 public:
  /** @brief The term coeff * lhs * rhs. */
  struct Product {
    double coeff;
    LinExpr lhs;
    LinExpr rhs;
  };

  QuadExpr(double constant = 0.0) : linear_(constant) {}  // NOLINT: implicit like GRBQuadExpr
  explicit QuadExpr(const LinExpr& linear) : linear_(linear) {}

  QuadExpr& operator+=(const QuadExpr& other);
  QuadExpr& operator+=(const LinExpr& other);
  QuadExpr& operator+=(double constant);
  QuadExpr& operator*=(double factor);

  /** @brief Add the term coeff * var1 * var2. */
  void addTerm(double coeff, int var1, int var2);
  /** @brief Add the term coeff * lhs * rhs. */
  void addProduct(double coeff, const LinExpr& lhs, const LinExpr& rhs);

  /** @brief Affine part, excluding the expansion of the products. */
  const LinExpr& linear() const { return linear_; }
  size_t numQuadTerms() const { return qvars1_.size(); }
  int getQuadVar1(size_t i) const { return qvars1_[i]; }
  int getQuadVar2(size_t i) const { return qvars2_[i]; }
  double getQuadCoeff(size_t i) const { return qcoeffs_[i]; }
  const std::vector<Product>& products() const { return products_; }

 private:
  LinExpr linear_;
  std::vector<int> qvars1_;
  std::vector<int> qvars2_;
  std::vector<double> qcoeffs_;
  std::vector<Product> products_;
};

QuadExpr operator*(const LinExpr& lhs, const LinExpr& rhs);
QuadExpr operator+(QuadExpr lhs, const QuadExpr& rhs);
QuadExpr operator+(QuadExpr lhs, const LinExpr& rhs);
QuadExpr operator+(QuadExpr lhs, double rhs);
QuadExpr operator*(double factor, QuadExpr expr);
QuadExpr operator*(QuadExpr expr, double factor);

/** @brief Linear constraint "expr sense 0", produced by comparing two LinExpr (cf. GRBTempConstr).
 */
struct TempConstr {
  LinExpr expr;
  char sense;
};

/** @brief Quadratic constraint "expr sense 0", produced by comparing a QuadExpr with a value. */
struct TempQConstr {
  QuadExpr expr;
  char sense;
};

TempConstr operator<=(const LinExpr& lhs, const LinExpr& rhs);
TempConstr operator>=(const LinExpr& lhs, const LinExpr& rhs);
TempConstr operator==(const LinExpr& lhs, const LinExpr& rhs);
TempQConstr operator<=(const QuadExpr& lhs, double rhs);
TempQConstr operator>=(const QuadExpr& lhs, double rhs);

/** @brief Handle to a linear constraint (cf. GRBConstr). */
struct Constr {
  int index{-1};
};
/** @brief Handle to a quadratic constraint (cf. GRBQConstr). */
struct QConstr {
  int index{-1};
};
/** @brief Handle to an indicator constraint (cf. GRBGenConstr). */
struct GenConstr {
  int index{-1};
};

/**
 * @brief Optimization model with Gurobi-like incremental editing, solved with HiGHS or SHOT.
 *
 * Variables and constraints live in slots that are recycled after removal, so repeatedly removing
 * and re-adding constraints between solves (as SolverGurobi does) does not grow memory. A Model
 * is not thread-safe, but separate Model instances can be solved concurrently from different
 * threads. HiGHS solves run in parallel; SHOT solves are serialized process-wide, since SHOT keeps
 * global state in its root search.
 */
class Model {
 public:
  Model() = default;

  /** @brief Add a variable with bounds [lb, ub], linear objective coefficient obj and type. */
  Var addVar(double lb, double ub, double obj, char type, const std::string& name = "");
  /** @brief Add a linear constraint built with <=, >= or == on LinExpr. */
  Constr addConstr(const TempConstr& constr, const std::string& name = "");
  /** @brief Add a quadratic constraint built with <= or >= on QuadExpr. */
  QConstr addQConstr(const TempQConstr& constr, const std::string& name = "");
  /** @brief Add the indicator constraint (binvar == binval) -> (expr sense rhs). */
  GenConstr addGenConstrIndicator(
      const Var& binvar, int binval, const LinExpr& expr, char sense, double rhs,
      const std::string& name = "");

  void remove(const Var& var);
  void remove(const Constr& constr);
  void remove(const QConstr& constr);
  void remove(const GenConstr& constr);

  /** @brief Set the objective function (replaces any previous objective). */
  void setObjective(const QuadExpr& objective, int sense = kMinimize);

  /** @brief Number of variables currently in the model. */
  int numVars() const;
  /** @brief Handles of all variables currently in the model. */
  std::vector<Var> getVars() const;

  // Parameters
  void setModelName(const std::string& name) { model_name_ = name; }
  void setOutputFlag(bool output) { output_flag_ = output; }
  /** @brief Wall-clock time limit for optimize() in seconds. */
  void setTimeLimit(double seconds) { time_limit_ = seconds; }
  /** @brief Threads used by HiGHS and by SHOT's MIP subsolver (0 = solver default). */
  void setThreads(int threads) { threads_ = threads; }
  /** @brief Relative and absolute optimality gaps for mixed-integer problems. */
  void setMipGap(double relative, double absolute) {
    mip_rel_gap_ = relative;
    mip_abs_gap_ = absolute;
  }
  /** @brief Flag polled during optimization; when true the solve is abandoned (cf. GRBCallback).
   */
  void setTerminationFlag(const std::atomic<bool>* flag) { terminate_flag_ = flag; }

  /** @brief Solve the model. Errors are reported through status() (kError), not exceptions. */
  Status optimize();

  Status status() const { return status_; }
  /** @brief Objective value of the last solution (NaN if none). */
  double objVal() const { return obj_val_; }
  /** @brief Wall-clock time of the last optimize() call in seconds. */
  double runtime() const { return runtime_; }
  Backend lastBackend() const { return last_backend_; }
  /** @brief Message describing the last error, if any. */
  const std::string& lastError() const { return last_error_; }

  /** @brief Value of a variable in the last solution (throws SolverException if none). */
  double value(const Var& var) const;
  /** @brief Value of a linear expression in the last solution. */
  double value(const LinExpr& expr) const;

 private:
  struct VarData {
    double lb;
    double ub;
    double obj;
    char type;
    std::string name;
    bool active;
  };
  struct LinConData {
    LinExpr expr;
    char sense;
    std::string name;
    bool active;
  };
  struct QuadConData {
    QuadExpr expr;
    char sense;
    std::string name;
    bool active;
  };
  struct IndicatorData {
    int binvar;
    int binval;
    LinExpr expr;  // expr - rhs, so the implied constraint is "expr sense 0"
    char sense;
    std::string name;
    bool active;
  };

  // A model assembled for the backends, in terms of the active (column-indexed) variables.
  struct Assembled;

  template <typename T>
  static int allocateSlot(std::vector<T>& slots, std::vector<int>& free_slots, T&& data);

  bool shouldTerminate() const { return terminate_flag_ && terminate_flag_->load(); }
  Status assemble(Assembled& a) const;
  Status solveWithHighs(const Assembled& a, double time_limit);
  Status solveWithShot(const Assembled& a, double time_limit);
  void polishSolution(const Assembled& a, const std::chrono::steady_clock::time_point& start);
  void storeSolution(const Assembled& a, const std::vector<double>& col_values);
  static double objectiveValue(const Assembled& a, const std::vector<double>& col_values);

  std::vector<VarData> vars_;
  std::vector<int> free_vars_;
  std::vector<LinConData> lin_cons_;
  std::vector<int> free_lin_cons_;
  std::vector<QuadConData> quad_cons_;
  std::vector<int> free_quad_cons_;
  std::vector<IndicatorData> indicators_;
  std::vector<int> free_indicators_;
  QuadExpr objective_;
  int objective_sense_{kMinimize};

  std::string model_name_;
  bool output_flag_{false};
  double time_limit_{kInf};
  int threads_{0};
  double mip_rel_gap_{1e-4};  // Gurobi's default MIPGap
  double mip_abs_gap_{1e-6};
  const std::atomic<bool>* terminate_flag_{nullptr};

  Status status_{Status::kLoaded};
  double obj_val_{std::numeric_limits<double>::quiet_NaN()};
  double runtime_{0.0};
  Backend last_backend_{Backend::kNone};
  std::string last_error_;
  std::vector<double> solution_;  // indexed by variable slot
};

}  // namespace hs
