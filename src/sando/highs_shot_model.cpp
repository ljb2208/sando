/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <sando/highs_shot_model.hpp>

#include <dlfcn.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <utility>

#include "highs_shot_assembled.hpp"

// HiGHS (SHOT is only used through the backend plugin, see highs_shot_assembled.hpp)
#include "Highs.h"

// Infinite bounds are represented by IEEE infinity, which -ffast-math assumes away
#ifdef __FAST_MATH__
#error "highs_shot_model.cpp must be compiled without -ffast-math (add -fno-fast-math)"
#endif

namespace hs {

namespace {

// Coefficients smaller than this are dropped when assembling rows. This matches HiGHS'
// default small_matrix_value, so both backends see the same model.
constexpr double kCoeffTol = 1e-9;
// Feasibility tolerance used to check rows left without variables (Gurobi's FeasibilityTol).
constexpr double kFeasTol = 1e-6;
// Big-M used for an indicator constraint whose expression is unbounded over the variable bounds
// even after bound tightening.
constexpr double kBigMFallback = 1e5;

using Clock = std::chrono::steady_clock;

double secondsSince(const Clock::time_point& start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

using detail::senseToBounds;
using detail::SparseAffine;
using detail::SquareTerm;

// HiGHS' global scheduler is thread-local and keeps the number of threads it was created with; a
// solve that asks for another (non-zero) number fails. Reset it when this thread last used a
// different number.
void prepareHighsScheduler(int threads) {
  thread_local int scheduler_threads = -1;
  if (threads <= 0) return;
  if (scheduler_threads != -1 && scheduler_threads != threads) Highs::resetGlobalScheduler(true);
  scheduler_threads = threads;
}

// SHOT solves are serialized process-wide: SHOT keeps the state of its root search in global
// variables (RootsearchMethodBoost.cpp), so concurrent solves corrupt each other.
std::timed_mutex& shotMutex() {
  static std::timed_mutex mutex;
  return mutex;
}

// Loads the SHOT backend plugin (once) and returns its entry point, or nullptr with `error` set.
// It is looked up at $SANDO_SHOT_BACKEND, then at its install location, then on the library path.
detail::ShotSolveFunction shotBackend(std::string& error) {
  static std::once_flag once;
  static detail::ShotSolveFunction solve = nullptr;
  static std::string load_error;
  std::call_once(once, [] {
    std::vector<std::string> candidates;
    if (const char* path = std::getenv("SANDO_SHOT_BACKEND")) candidates.emplace_back(path);
#ifdef SANDO_SHOT_BACKEND_PATH
    candidates.emplace_back(SANDO_SHOT_BACKEND_PATH);
#endif
    candidates.emplace_back("libsando_shot_backend.so");
    for (const auto& candidate : candidates) {
      // RTLD_LOCAL keeps SHOT's symbols out of the global scope, RTLD_DEEPBIND makes SHOT bind to
      // its own symbols first (see highs_shot_assembled.hpp)
      void* handle = dlopen(candidate.c_str(), RTLD_NOW | RTLD_LOCAL | RTLD_DEEPBIND);
      if (!handle) {
        load_error += std::string(dlerror()) + "; ";
        continue;
      }
      solve = reinterpret_cast<detail::ShotSolveFunction>(dlsym(handle, detail::kShotSolveSymbol));
      if (solve) return;  // the plugin stays loaded for the lifetime of the process
      load_error += candidate + ": no symbol " + detail::kShotSolveSymbol + "; ";
      dlclose(handle);
    }
  });
  if (!solve) error = "could not load the SHOT backend plugin: " + load_error;
  return solve;
}

}  // namespace

// --------------------------------------------------------------------------------------------
// LinExpr
// --------------------------------------------------------------------------------------------

LinExpr::LinExpr(const Var& var, double coeff) {
  if (!var.valid()) throw SolverException("LinExpr: invalid variable");
  vars_.push_back(var.index());
  coeffs_.push_back(coeff);
}

LinExpr& LinExpr::operator+=(const LinExpr& other) {
  if (&other == this) return *this *= 2.0;
  vars_.insert(vars_.end(), other.vars_.begin(), other.vars_.end());
  coeffs_.insert(coeffs_.end(), other.coeffs_.begin(), other.coeffs_.end());
  constant_ += other.constant_;
  return *this;
}

LinExpr& LinExpr::operator-=(const LinExpr& other) {
  if (&other == this) return *this *= 0.0;
  vars_.insert(vars_.end(), other.vars_.begin(), other.vars_.end());
  for (double c : other.coeffs_) coeffs_.push_back(-c);
  constant_ -= other.constant_;
  return *this;
}

LinExpr& LinExpr::operator*=(double factor) {
  for (double& c : coeffs_) c *= factor;
  constant_ *= factor;
  return *this;
}

LinExpr& LinExpr::operator/=(double divisor) {
  for (double& c : coeffs_) c /= divisor;
  constant_ /= divisor;
  return *this;
}

LinExpr operator+(LinExpr lhs, const LinExpr& rhs) { return lhs += rhs; }
LinExpr operator-(LinExpr lhs, const LinExpr& rhs) { return lhs -= rhs; }
LinExpr operator+(LinExpr lhs, double rhs) { return lhs += LinExpr(rhs); }
LinExpr operator+(double lhs, LinExpr rhs) { return rhs += LinExpr(lhs); }
LinExpr operator-(LinExpr lhs, double rhs) { return lhs -= LinExpr(rhs); }
LinExpr operator-(double lhs, LinExpr rhs) { return (rhs *= -1.0) += LinExpr(lhs); }
LinExpr operator-(LinExpr expr) { return expr *= -1.0; }
LinExpr operator*(LinExpr expr, double factor) { return expr *= factor; }
LinExpr operator*(double factor, LinExpr expr) { return expr *= factor; }
LinExpr operator/(LinExpr expr, double divisor) { return expr /= divisor; }

// --------------------------------------------------------------------------------------------
// QuadExpr
// --------------------------------------------------------------------------------------------

QuadExpr& QuadExpr::operator+=(const QuadExpr& other) {
  if (&other == this) return *this *= 2.0;
  linear_ += other.linear_;
  qvars1_.insert(qvars1_.end(), other.qvars1_.begin(), other.qvars1_.end());
  qvars2_.insert(qvars2_.end(), other.qvars2_.begin(), other.qvars2_.end());
  qcoeffs_.insert(qcoeffs_.end(), other.qcoeffs_.begin(), other.qcoeffs_.end());
  products_.insert(products_.end(), other.products_.begin(), other.products_.end());
  return *this;
}

QuadExpr& QuadExpr::operator+=(const LinExpr& other) {
  linear_ += other;
  return *this;
}

QuadExpr& QuadExpr::operator+=(double constant) {
  linear_ += LinExpr(constant);
  return *this;
}

QuadExpr& QuadExpr::operator*=(double factor) {
  linear_ *= factor;
  for (double& c : qcoeffs_) c *= factor;
  for (auto& p : products_) p.coeff *= factor;
  return *this;
}

void QuadExpr::addTerm(double coeff, int var1, int var2) {
  qvars1_.push_back(var1);
  qvars2_.push_back(var2);
  qcoeffs_.push_back(coeff);
}

void QuadExpr::addProduct(double coeff, const LinExpr& lhs, const LinExpr& rhs) {
  products_.push_back({coeff, lhs, rhs});
}

QuadExpr operator*(const LinExpr& lhs, const LinExpr& rhs) {
  QuadExpr result;
  result.addProduct(1.0, lhs, rhs);
  return result;
}

QuadExpr operator+(QuadExpr lhs, const QuadExpr& rhs) { return lhs += rhs; }
QuadExpr operator+(QuadExpr lhs, const LinExpr& rhs) { return lhs += rhs; }
QuadExpr operator+(QuadExpr lhs, double rhs) { return lhs += rhs; }
QuadExpr operator*(double factor, QuadExpr expr) { return expr *= factor; }
QuadExpr operator*(QuadExpr expr, double factor) { return expr *= factor; }

// --------------------------------------------------------------------------------------------
// Constraint builders
// --------------------------------------------------------------------------------------------

TempConstr operator<=(const LinExpr& lhs, const LinExpr& rhs) { return {lhs - rhs, kLessEqual}; }
TempConstr operator>=(const LinExpr& lhs, const LinExpr& rhs) {
  return {lhs - rhs, kGreaterEqual};
}
TempConstr operator==(const LinExpr& lhs, const LinExpr& rhs) { return {lhs - rhs, kEqual}; }
TempQConstr operator<=(const QuadExpr& lhs, double rhs) { return {lhs + (-rhs), kLessEqual}; }
TempQConstr operator>=(const QuadExpr& lhs, double rhs) { return {lhs + (-rhs), kGreaterEqual}; }

// --------------------------------------------------------------------------------------------
// Model: editing
// --------------------------------------------------------------------------------------------

template <typename T>
int Model::allocateSlot(std::vector<T>& slots, std::vector<int>& free_slots, T&& data) {
  if (!free_slots.empty()) {
    const int index = free_slots.back();
    free_slots.pop_back();
    slots[index] = std::move(data);
    return index;
  }
  slots.push_back(std::move(data));
  return static_cast<int>(slots.size()) - 1;
}

Var Model::addVar(double lb, double ub, double obj, char type, const std::string& name) {
  if (type != kContinuous && type != kBinary)
    throw SolverException("Model::addVar: unsupported variable type");
  return Var(allocateSlot(vars_, free_vars_, VarData{lb, ub, obj, type, name, true}));
}

Constr Model::addConstr(const TempConstr& constr, const std::string& name) {
  return Constr{
      allocateSlot(lin_cons_, free_lin_cons_, LinConData{constr.expr, constr.sense, name, true})};
}

QConstr Model::addQConstr(const TempQConstr& constr, const std::string& name) {
  return QConstr{
      allocateSlot(quad_cons_, free_quad_cons_, QuadConData{constr.expr, constr.sense, name, true})};
}

GenConstr Model::addGenConstrIndicator(
    const Var& binvar, int binval, const LinExpr& expr, char sense, double rhs,
    const std::string& name) {
  if (!binvar.valid() || binvar.index() >= static_cast<int>(vars_.size()) ||
      !vars_[binvar.index()].active || vars_[binvar.index()].type != kBinary)
    throw SolverException("Model::addGenConstrIndicator: indicator variable must be binary");
  if (binval != 0 && binval != 1)
    throw SolverException("Model::addGenConstrIndicator: indicator value must be 0 or 1");
  return GenConstr{allocateSlot(
      indicators_, free_indicators_,
      IndicatorData{binvar.index(), binval, expr - LinExpr(rhs), sense, name, true})};
}

void Model::remove(const Var& var) {
  const int index = var.index();
  if (index < 0 || index >= static_cast<int>(vars_.size()) || !vars_[index].active) return;

  // As in Gurobi, removing a variable removes its terms from every constraint and the objective,
  // and removes the indicator constraints it controls. This also prevents a recycled slot from
  // being referenced by a stale term.
  auto purge = [index](const LinExpr& expr) {
    LinExpr result(expr.getConstant());
    for (size_t k = 0; k < expr.size(); ++k)
      if (expr.getVarIndex(k) != index) result += LinExpr(Var(expr.getVarIndex(k)), expr.getCoeff(k));
    return result;
  };
  auto purgeQuad = [&purge, index](const QuadExpr& expr) {
    QuadExpr result(purge(expr.linear()));
    for (size_t k = 0; k < expr.numQuadTerms(); ++k)
      if (expr.getQuadVar1(k) != index && expr.getQuadVar2(k) != index)
        result.addTerm(expr.getQuadCoeff(k), expr.getQuadVar1(k), expr.getQuadVar2(k));
    for (const auto& p : expr.products()) result.addProduct(p.coeff, purge(p.lhs), purge(p.rhs));
    return result;
  };
  auto references = [index](const LinExpr& expr) {
    for (size_t k = 0; k < expr.size(); ++k)
      if (expr.getVarIndex(k) == index) return true;
    return false;
  };
  auto referencesQuad = [&references, index](const QuadExpr& expr) {
    if (references(expr.linear())) return true;
    for (size_t k = 0; k < expr.numQuadTerms(); ++k)
      if (expr.getQuadVar1(k) == index || expr.getQuadVar2(k) == index) return true;
    for (const auto& p : expr.products())
      if (references(p.lhs) || references(p.rhs)) return true;
    return false;
  };

  for (auto& c : lin_cons_)
    if (c.active && references(c.expr)) c.expr = purge(c.expr);
  for (auto& c : quad_cons_)
    if (c.active && referencesQuad(c.expr)) c.expr = purgeQuad(c.expr);
  for (size_t i = 0; i < indicators_.size(); ++i) {
    auto& c = indicators_[i];
    if (!c.active) continue;
    if (c.binvar == index)
      remove(GenConstr{static_cast<int>(i)});
    else if (references(c.expr))
      c.expr = purge(c.expr);
  }
  if (referencesQuad(objective_)) objective_ = purgeQuad(objective_);

  vars_[index].active = false;
  vars_[index].name.clear();
  free_vars_.push_back(index);
}

void Model::remove(const Constr& constr) {
  const int i = constr.index;
  if (i < 0 || i >= static_cast<int>(lin_cons_.size()) || !lin_cons_[i].active) return;
  lin_cons_[i] = LinConData{LinExpr(), kEqual, std::string(), false};
  free_lin_cons_.push_back(i);
}

void Model::remove(const QConstr& constr) {
  const int i = constr.index;
  if (i < 0 || i >= static_cast<int>(quad_cons_.size()) || !quad_cons_[i].active) return;
  quad_cons_[i] = QuadConData{QuadExpr(), kEqual, std::string(), false};
  free_quad_cons_.push_back(i);
}

void Model::remove(const GenConstr& constr) {
  const int i = constr.index;
  if (i < 0 || i >= static_cast<int>(indicators_.size()) || !indicators_[i].active) return;
  indicators_[i] = IndicatorData{-1, 1, LinExpr(), kEqual, std::string(), false};
  free_indicators_.push_back(i);
}

void Model::setObjective(const QuadExpr& objective, int sense) {
  objective_ = objective;
  objective_sense_ = sense;
}

int Model::numVars() const { return static_cast<int>(vars_.size() - free_vars_.size()); }

std::vector<Var> Model::getVars() const {
  std::vector<Var> result;
  result.reserve(numVars());
  for (size_t i = 0; i < vars_.size(); ++i)
    if (vars_[i].active) result.emplace_back(static_cast<int>(i));
  return result;
}

double Model::value(const Var& var) const {
  const int i = var.index();
  if (status_ != Status::kOptimal || i < 0 || i >= static_cast<int>(solution_.size()) ||
      std::isnan(solution_[i]))
    throw SolverException("Model::value: no solution value available for this variable");
  return solution_[i];
}

double Model::value(const LinExpr& expr) const {
  double result = expr.getConstant();
  for (size_t k = 0; k < expr.size(); ++k)
    result += expr.getCoeff(k) * value(Var(expr.getVarIndex(k)));
  return result;
}

// --------------------------------------------------------------------------------------------
// Model: assembly
// --------------------------------------------------------------------------------------------


namespace {

// Map a LinExpr to columns, merging repeated variables and dropping negligible coefficients.
// `dense` must have one zero-initialized entry per column and is left zeroed.
bool compress(
    const LinExpr& expr, const std::vector<int>& slot_col, std::vector<double>& dense,
    SparseAffine& out) {
  out.cols.clear();
  out.vals.clear();
  out.constant = expr.getConstant();
  for (size_t k = 0; k < expr.size(); ++k) {
    const int slot = expr.getVarIndex(k);
    if (slot < 0 || slot >= static_cast<int>(slot_col.size()) || slot_col[slot] < 0) return false;
    const int col = slot_col[slot];
    if (dense[col] == 0.0) out.cols.push_back(col);
    dense[col] += expr.getCoeff(k);
    if (dense[col] == 0.0) dense[col] = 1e-300;  // keep the column registered once
  }
  size_t kept = 0;
  for (size_t k = 0; k < out.cols.size(); ++k) {
    const int col = out.cols[k];
    const double v = dense[col];
    dense[col] = 0.0;
    if (std::abs(v) > kCoeffTol) {
      out.cols[kept++] = col;
      out.vals.push_back(v);
    }
  }
  out.cols.resize(kept);
  return true;
}

// Quadratic terms coeff * x_var1 * x_var2 (by variable slot)
struct QuadTerms {
  std::vector<int> var1, var2;
  std::vector<double> coeff;
  void add(double c, int v1, int v2) {
    var1.push_back(v1);
    var2.push_back(v2);
    coeff.push_back(c);
  }
};

// Expand the products of a quadratic expression into its full affine part and quadratic terms.
void expandQuad(const QuadExpr& expr, LinExpr& linear, QuadTerms& quad) {
  linear = expr.linear();
  quad = QuadTerms();
  for (size_t k = 0; k < expr.numQuadTerms(); ++k)
    quad.add(expr.getQuadCoeff(k), expr.getQuadVar1(k), expr.getQuadVar2(k));
  // c * (sum a_i x_i + p)(sum b_j x_j + q) = c sum a_i b_j x_i x_j + c q sum a_i x_i
  //                                          + c p sum b_j x_j + c p q
  for (const auto& prod : expr.products()) {
    const LinExpr& l = prod.lhs;
    const LinExpr& r = prod.rhs;
    for (size_t i = 0; i < l.size(); ++i)
      for (size_t j = 0; j < r.size(); ++j)
        quad.add(prod.coeff * l.getCoeff(i) * r.getCoeff(j), l.getVarIndex(i), r.getVarIndex(j));
    for (size_t i = 0; i < l.size(); ++i)
      linear += LinExpr(Var(l.getVarIndex(i)), prod.coeff * l.getCoeff(i) * r.getConstant());
    for (size_t j = 0; j < r.size(); ++j)
      linear += LinExpr(Var(r.getVarIndex(j)), prod.coeff * r.getCoeff(j) * l.getConstant());
    linear += LinExpr(prod.coeff * l.getConstant() * r.getConstant());
  }
}

bool sameExpr(const LinExpr& a, const LinExpr& b) {
  if (a.size() != b.size() || a.getConstant() != b.getConstant()) return false;
  for (size_t k = 0; k < a.size(); ++k)
    if (a.getVarIndex(k) != b.getVarIndex(k) || a.getCoeff(k) != b.getCoeff(k)) return false;
  return true;
}

// If the quadratic part of `expr` is a nonnegative combination of squares of affine expressions,
// return them; the remaining affine part is then expr.linear().
bool extractSquares(
    const QuadExpr& expr, const std::vector<int>& slot_col, std::vector<double>& dense,
    std::vector<SquareTerm>& squares) {
  squares.clear();
  if (expr.numQuadTerms() > 0) return false;
  for (const auto& prod : expr.products()) {
    if (prod.coeff < 0.0 || !sameExpr(prod.lhs, prod.rhs)) return false;
    if (prod.coeff == 0.0) continue;
    SquareTerm square{prod.coeff, {}};
    if (!compress(prod.lhs, slot_col, dense, square.expr)) return false;
    squares.push_back(std::move(square));
  }
  return !squares.empty();
}

// Merge quadratic terms into (i <= j) column pairs, dropping negligible coefficients.
bool compressQuad(
    const QuadTerms& quad, const std::vector<int>& slot_col, std::vector<int>& qi,
    std::vector<int>& qj, std::vector<double>& qv) {
  std::vector<std::pair<std::pair<int, int>, double>> terms;
  terms.reserve(quad.coeff.size());
  for (size_t k = 0; k < quad.coeff.size(); ++k) {
    const int s1 = quad.var1[k], s2 = quad.var2[k];
    if (s1 < 0 || s2 < 0 || s1 >= static_cast<int>(slot_col.size()) ||
        s2 >= static_cast<int>(slot_col.size()) || slot_col[s1] < 0 || slot_col[s2] < 0)
      return false;
    const int c1 = std::min(slot_col[s1], slot_col[s2]);
    const int c2 = std::max(slot_col[s1], slot_col[s2]);
    terms.push_back({{c1, c2}, quad.coeff[k]});
  }
  std::sort(terms.begin(), terms.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  qi.clear();
  qj.clear();
  qv.clear();
  for (size_t k = 0; k < terms.size();) {
    const auto key = terms[k].first;
    double sum = 0.0;
    for (; k < terms.size() && terms[k].first == key; ++k) sum += terms[k].second;
    if (std::abs(sum) > kCoeffTol) {
      qi.push_back(key.first);
      qj.push_back(key.second);
      qv.push_back(sum);
    }
  }
  return true;
}

// Feasibility-based bound tightening over the linear rows. Only used to derive big-M constants,
// so the bounds are relaxed slightly to stay valid under round-off.
void tightenBounds(
    const std::vector<int>& row_start, const std::vector<int>& row_index,
    const std::vector<double>& row_value, const std::vector<double>& row_lo,
    const std::vector<double>& row_up, std::vector<double>& lb, std::vector<double>& ub) {
  const int num_rows = static_cast<int>(row_lo.size());
  const int max_passes = 20;
  for (int pass = 0; pass < max_passes; ++pass) {
    bool changed = false;
    for (int r = 0; r < num_rows; ++r) {
      // Activity bounds of the row, counting infinite contributions separately
      double min_act = 0.0, max_act = 0.0;
      int min_inf = 0, max_inf = 0;
      for (int k = row_start[r]; k < row_start[r + 1]; ++k) {
        const double a = row_value[k];
        const int j = row_index[k];
        const double cmin = a > 0 ? a * lb[j] : a * ub[j];
        const double cmax = a > 0 ? a * ub[j] : a * lb[j];
        if (std::isinf(cmin)) ++min_inf; else min_act += cmin;
        if (std::isinf(cmax)) ++max_inf; else max_act += cmax;
      }
      if (min_inf > 1 && max_inf > 1) continue;

      for (int k = row_start[r]; k < row_start[r + 1]; ++k) {
        const double a = row_value[k];
        const int j = row_index[k];
        if (std::abs(a) < 1e-7) continue;
        const double cmin = a > 0 ? a * lb[j] : a * ub[j];
        const double cmax = a > 0 ? a * ub[j] : a * lb[j];

        // Minimum / maximum activity of the other terms of the row
        double others_min = -kInf, others_max = kInf;
        if (std::isinf(cmin)) {
          if (min_inf == 1) others_min = min_act;
        } else if (min_inf == 0) {
          others_min = min_act - cmin;
        }
        if (std::isinf(cmax)) {
          if (max_inf == 1) others_max = max_act;
        } else if (max_inf == 0) {
          others_max = max_act - cmax;
        }

        double new_lb = -kInf, new_ub = kInf;
        if (!std::isinf(row_up[r]) && !std::isinf(others_min)) {
          const double bound = (row_up[r] - others_min) / a;  // a * x_j <= row_up - others_min
          if (a > 0) new_ub = bound; else new_lb = bound;
        }
        if (!std::isinf(row_lo[r]) && !std::isinf(others_max)) {
          const double bound = (row_lo[r] - others_max) / a;  // a * x_j >= row_lo - others_max
          if (a > 0) new_lb = std::max(new_lb, bound); else new_ub = std::min(new_ub, bound);
        }

        if (!std::isinf(new_ub)) {
          new_ub += 1e-7 * (1.0 + std::abs(new_ub));
          if (new_ub < ub[j] - 1e-6 * (1.0 + std::abs(ub[j])) || std::isinf(ub[j])) {
            ub[j] = new_ub;
            changed = true;
          }
        }
        if (!std::isinf(new_lb)) {
          new_lb -= 1e-7 * (1.0 + std::abs(new_lb));
          if (new_lb > lb[j] + 1e-6 * (1.0 + std::abs(lb[j])) || std::isinf(lb[j])) {
            lb[j] = new_lb;
            changed = true;
          }
        }
      }
    }
    if (!changed) break;
  }
}

}  // namespace

Status Model::assemble(Assembled& a) const {
  // Columns
  a.slot_col.assign(vars_.size(), -1);
  for (size_t s = 0; s < vars_.size(); ++s) {
    const VarData& v = vars_[s];
    if (!v.active) continue;
    a.slot_col[s] = a.num_cols++;
    a.col_slot.push_back(static_cast<int>(s));
    double lb = v.lb, ub = v.ub;
    if (v.type == kBinary) {
      lb = std::ceil(std::max(lb, 0.0) - kFeasTol);
      ub = std::floor(std::min(ub, 1.0) + kFeasTol);
    }
    a.lb.push_back(lb);
    a.ub.push_back(ub);
    a.cost.push_back(v.obj);
    a.is_binary.push_back(v.type == kBinary);
  }

  std::vector<double> dense(a.num_cols, 0.0);
  SparseAffine aff;

  // Objective
  a.sense = objective_sense_;
  LinExpr linear;
  QuadTerms quad;
  expandQuad(objective_, linear, quad);
  if (!compress(linear, a.slot_col, dense, aff) ||
      !compressQuad(quad, a.slot_col, a.obj_qi, a.obj_qj, a.obj_qv))
    throw SolverException("objective references a removed variable");
  a.obj_sos_cost = a.cost;
  for (size_t k = 0; k < aff.cols.size(); ++k) a.cost[aff.cols[k]] += aff.vals[k];
  a.offset = aff.constant;
  a.obj_is_sos = extractSquares(objective_, a.slot_col, dense, a.obj_squares);
  if (a.obj_is_sos) {
    if (!compress(objective_.linear(), a.slot_col, dense, aff))
      throw SolverException("objective references a removed variable");
    for (size_t k = 0; k < aff.cols.size(); ++k) a.obj_sos_cost[aff.cols[k]] += aff.vals[k];
    a.obj_sos_offset = aff.constant;
  }

  // Adds "aff sense 0" as a linear row; returns false if it is a violated constant row.
  auto addLinear = [&a](const SparseAffine& row, char sense) {
    double lo, up;
    senseToBounds(sense, row.constant, lo, up);
    if (row.cols.empty()) return lo <= kFeasTol && up >= -kFeasTol;
    a.addRow(row.cols, row.vals, lo, up);
    return true;
  };

  bool infeasible = false;

  // Linear constraints
  for (const auto& c : lin_cons_) {
    if (!c.active) continue;
    if (!compress(c.expr, a.slot_col, dense, aff))
      throw SolverException("constraint " + c.name + " references a removed variable");
    infeasible |= !addLinear(aff, c.sense);
  }

  // Quadratic constraints (those whose quadratic part cancels are linear)
  for (const auto& c : quad_cons_) {
    if (!c.active) continue;
    Assembled::QRow q;
    q.sense = c.sense;
    expandQuad(c.expr, linear, quad);
    if (!compress(linear, a.slot_col, dense, q.lin) ||
        !compressQuad(quad, a.slot_col, q.qi, q.qj, q.qv))
      throw SolverException("constraint " + c.name + " references a removed variable");
    if (q.qv.empty()) {
      infeasible |= !addLinear(q.lin, c.sense);
      continue;
    }
    q.is_sos = extractSquares(c.expr, a.slot_col, dense, q.squares) &&
               compress(c.expr.linear(), a.slot_col, dense, q.sos_lin);
    a.qrows.push_back(std::move(q));
  }

  // Binary variables fixed by single-variable rows, e.g. "sum of one binary >= 1" or "b == 0"
  for (int r = 0; r < a.numRows(); ++r) {
    if (a.row_start[r + 1] - a.row_start[r] != 1) continue;
    const int j = a.row_index[a.row_start[r]];
    const double v = a.row_value[a.row_start[r]];
    if (!a.is_binary[j]) continue;
    double lo = a.row_lo[r] / v, up = a.row_up[r] / v;
    if (v < 0) std::swap(lo, up);
    if (!std::isinf(lo)) a.lb[j] = std::max(a.lb[j], std::ceil(lo - kFeasTol));
    if (!std::isinf(up)) a.ub[j] = std::min(a.ub[j], std::floor(up + kFeasTol));
  }
  for (int j = 0; j < a.num_cols; ++j) {
    if (a.lb[j] > a.ub[j] + kFeasTol) infeasible = true;
    if (a.is_binary[j] && a.lb[j] < a.ub[j]) a.has_integer = true;
  }

  // Indicator constraints: enforced or dropped when their binary is fixed, big-M otherwise
  struct PendingIndicator {
    int bincol;
    int binval;
    SparseAffine expr;
    char sense;
  };
  std::vector<PendingIndicator> pending;
  for (const auto& c : indicators_) {
    if (!c.active) continue;
    if (c.binvar < 0 || c.binvar >= static_cast<int>(a.slot_col.size()) || a.slot_col[c.binvar] < 0 ||
        !compress(c.expr, a.slot_col, dense, aff))
      throw SolverException("indicator constraint " + c.name + " references a removed variable");
    const int bincol = a.slot_col[c.binvar];
    if (a.lb[bincol] == a.ub[bincol]) {
      if (static_cast<int>(a.lb[bincol]) == c.binval) infeasible |= !addLinear(aff, c.sense);
      continue;
    }
    pending.push_back({bincol, c.binval, aff, c.sense});
  }

  if (!pending.empty()) {
    std::vector<double> lb = a.lb, ub = a.ub;
    tightenBounds(a.row_start, a.row_index, a.row_value, a.row_lo, a.row_up, lb, ub);

    // Adds "expr <= 0 if bin == binval" as expr <= M * (1 - bin) (binval 1) or expr <= M * bin
    auto addBigM = [&](const PendingIndicator& ind, const SparseAffine& expr) {
      double max_act = expr.constant;
      for (size_t k = 0; k < expr.cols.size(); ++k) {
        const double v = expr.vals[k];
        max_act += v > 0 ? v * ub[expr.cols[k]] : v * lb[expr.cols[k]];
      }
      if (max_act <= 0.0) return;  // implied by the bounds
      double big_m = kBigMFallback;
      if (!std::isinf(max_act) && !std::isnan(max_act)) big_m = max_act * (1.0 + 1e-6) + 1e-6;
      std::vector<int> cols = expr.cols;
      std::vector<double> vals = expr.vals;
      cols.push_back(ind.bincol);
      if (ind.binval == 1) {
        vals.push_back(big_m);
        a.addRow(cols, vals, -kInf, big_m - expr.constant);
      } else {
        vals.push_back(-big_m);
        a.addRow(cols, vals, -kInf, -expr.constant);
      }
    };
    for (const auto& ind : pending) {
      SparseAffine negated = ind.expr;
      for (double& v : negated.vals) v = -v;
      negated.constant = -negated.constant;
      if (ind.sense == kLessEqual || ind.sense == kEqual) addBigM(ind, ind.expr);
      if (ind.sense == kGreaterEqual || ind.sense == kEqual) addBigM(ind, negated);
    }
  }

  return infeasible ? Status::kInfeasible : Status::kLoaded;
}

void Model::storeSolution(const Assembled& a, const std::vector<double>& col_values) {
  solution_.assign(vars_.size(), std::numeric_limits<double>::quiet_NaN());
  for (int j = 0; j < a.num_cols && j < static_cast<int>(col_values.size()); ++j)
    solution_[a.col_slot[j]] = col_values[j];
}

// --------------------------------------------------------------------------------------------
// Model: solving
// --------------------------------------------------------------------------------------------

Status Model::optimize() {
  const auto start = Clock::now();
  status_ = Status::kLoaded;
  obj_val_ = std::numeric_limits<double>::quiet_NaN();
  last_backend_ = Backend::kNone;
  last_error_.clear();
  solution_.clear();

  try {
    Assembled a;
    status_ = assemble(a);
    if (status_ == Status::kLoaded) {
      if (shouldTerminate()) {
        status_ = Status::kInterrupted;
      } else {
        const double remaining = time_limit_ - secondsSince(start);
        if (remaining <= 0.0) {
          status_ = Status::kTimeLimit;
        } else {
          // HiGHS handles LP, QP and MILP; SHOT handles MIQP and quadratic constraints.
          const bool quadratic = !a.obj_qv.empty() || !a.qrows.empty();
          const bool use_shot = !a.qrows.empty() || (a.has_integer && quadratic);
          status_ = use_shot ? solveWithShot(a, remaining) : solveWithHighs(a, remaining);
          if (use_shot && status_ == Status::kOptimal && a.has_integer) polishSolution(a, start);
        }
      }
    }
  } catch (const std::exception& e) {
    status_ = Status::kError;
    last_error_ = e.what();
  }

  if (status_ != Status::kOptimal) solution_.clear();
  runtime_ = secondsSince(start);
  return status_;
}

void Model::polishSolution(const Assembled& a, const std::chrono::steady_clock::time_point& start) {
  // SHOT's mixed-integer solution is accurate only to its integrality tolerance, which the big-M
  // rows of the indicator constraints amplify into constraint violations. With the binaries fixed
  // at SHOT's choice the big-M rows become exact, and the remaining continuous problem is solved
  // again: a convex QP exactly with HiGHS, or a convex QCQP with SHOT (which then uses Ipopt).
  const double remaining = time_limit_ - secondsSince(start);
  if (remaining <= 0.0 || shouldTerminate()) return;

  Assembled fixed = a;
  for (int j = 0; j < fixed.num_cols; ++j) {
    if (!fixed.is_binary[j]) continue;
    fixed.lb[j] = fixed.ub[j] = std::round(solution_[a.col_slot[j]]);
  }
  fixed.has_integer = false;

  const std::vector<double> shot_solution = solution_;
  const double shot_obj = obj_val_;
  const Status status =
      fixed.qrows.empty() ? solveWithHighs(fixed, remaining) : solveWithShot(fixed, remaining);
  if (status != Status::kOptimal) {
    solution_ = shot_solution;
    obj_val_ = shot_obj;
  }
  last_backend_ = Backend::kShot;
}

Status Model::solveWithHighs(const Assembled& a, double time_limit) {
  last_backend_ = Backend::kHighs;

  HighsModel model;
  HighsLp& lp = model.lp_;
  lp.model_name_ = model_name_;
  lp.num_col_ = a.num_cols;
  lp.num_row_ = a.numRows();
  lp.col_cost_ = a.cost;
  lp.col_lower_ = a.lb;
  lp.col_upper_ = a.ub;
  lp.row_lower_ = a.row_lo;
  lp.row_upper_ = a.row_up;
  lp.offset_ = a.offset;
  lp.sense_ = a.sense == kMaximize ? ObjSense::kMaximize : ObjSense::kMinimize;
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.num_col_ = a.num_cols;
  lp.a_matrix_.num_row_ = a.numRows();
  lp.a_matrix_.start_ = a.row_start;
  lp.a_matrix_.index_ = a.row_index;
  lp.a_matrix_.value_ = a.row_value;
  if (a.has_integer) {
    lp.integrality_.resize(a.num_cols, HighsVarType::kContinuous);
    for (int j = 0; j < a.num_cols; ++j)
      if (a.is_binary[j]) lp.integrality_[j] = HighsVarType::kInteger;
  }

  if (!a.obj_qv.empty()) {
    // HiGHS minimizes c'x + 0.5 x'Qx with the lower triangle of Q stored column-wise.
    // A term q * x_i * x_j (i < j) gives Q_ji = q; a term q * x_i^2 gives Q_ii = 2q.
    std::vector<std::vector<std::pair<int, double>>> columns(a.num_cols);
    for (size_t k = 0; k < a.obj_qv.size(); ++k) {
      const int i = a.obj_qi[k], j = a.obj_qj[k];
      columns[i].push_back({j, i == j ? 2.0 * a.obj_qv[k] : a.obj_qv[k]});
    }
    HighsHessian& hessian = model.hessian_;
    hessian.dim_ = a.num_cols;
    hessian.format_ = HessianFormat::kTriangular;
    hessian.start_.assign(1, 0);
    for (auto& column : columns) {
      std::sort(column.begin(), column.end());
      for (const auto& entry : column) {
        hessian.index_.push_back(entry.first);
        hessian.value_.push_back(entry.second);
      }
      hessian.start_.push_back(static_cast<HighsInt>(hessian.index_.size()));
    }
  }

  Highs highs;
  highs.setOptionValue("output_flag", output_flag_);
  highs.setOptionValue("threads", threads_);
  if (!std::isinf(time_limit)) highs.setOptionValue("time_limit", time_limit);
  if (a.has_integer) {
    highs.setOptionValue("mip_rel_gap", mip_rel_gap_);
    highs.setOptionValue("mip_abs_gap", mip_abs_gap_);
  }

  const std::atomic<bool>* flag = terminate_flag_;
  if (flag) {
    highs.setCallback(
        [flag](
            int, const std::string&, const HighsCallbackOutput*, HighsCallbackInput* data_in,
            void*) {
          if (data_in && flag->load()) data_in->user_interrupt = true;
        },
        nullptr);
    highs.startCallback(kCallbackSimplexInterrupt);
    highs.startCallback(kCallbackIpmInterrupt);
    highs.startCallback(kCallbackMipInterrupt);
  }

  if (highs.passModel(std::move(model)) == HighsStatus::kError) {
    last_error_ = "HiGHS rejected the model";
    return Status::kError;
  }

  prepareHighsScheduler(threads_);
  HighsStatus run_status = highs.run();
  if (run_status == HighsStatus::kError && highs.getModelStatus() == HighsModelStatus::kNotset) {
    // The thread's scheduler may have been created with another number of threads elsewhere
    Highs::resetGlobalScheduler(true);
    run_status = highs.run();
  }
  if (run_status == HighsStatus::kError) {
    last_error_ = "HiGHS run failed: " + highs.modelStatusToString(highs.getModelStatus());
    return Status::kError;
  }

  switch (highs.getModelStatus()) {
    case HighsModelStatus::kOptimal:
      storeSolution(a, highs.getSolution().col_value);
      obj_val_ = highs.getInfo().objective_function_value;
      return Status::kOptimal;
    case HighsModelStatus::kModelEmpty:
      storeSolution(a, std::vector<double>(a.num_cols, 0.0));
      obj_val_ = a.offset;
      return Status::kOptimal;
    case HighsModelStatus::kInfeasible:
      return Status::kInfeasible;
    case HighsModelStatus::kUnboundedOrInfeasible:
      return Status::kInfOrUnbd;
    case HighsModelStatus::kUnbounded:
      return Status::kUnbounded;
    case HighsModelStatus::kTimeLimit:
      return Status::kTimeLimit;
    case HighsModelStatus::kInterrupt:
      return Status::kInterrupted;
    case HighsModelStatus::kIterationLimit:
    case HighsModelStatus::kUnknown:
      return Status::kNumeric;
    default:
      last_error_ = "HiGHS: " + highs.modelStatusToString(highs.getModelStatus());
      return Status::kError;
  }
}

// Objective of the assembled (expanded) model at the given column values
double Model::objectiveValue(const Assembled& a, const std::vector<double>& x) {
  double value = a.offset;
  for (int j = 0; j < a.num_cols; ++j) value += a.cost[j] * x[j];
  for (size_t k = 0; k < a.obj_qv.size(); ++k) value += a.obj_qv[k] * x[a.obj_qi[k]] * x[a.obj_qj[k]];
  return value;
}

Status Model::solveWithShot(const Assembled& a, double time_limit) {
  last_backend_ = Backend::kShot;
  const auto start = Clock::now();

  // Wait for the process-wide SHOT lock, giving up when asked to stop or out of time
  std::unique_lock<std::timed_mutex> lock(shotMutex(), std::defer_lock);
  while (!lock.try_lock_for(std::chrono::milliseconds(1))) {
    if (shouldTerminate()) return Status::kInterrupted;
    if (secondsSince(start) >= time_limit) return Status::kTimeLimit;
  }
  if (shouldTerminate()) return Status::kInterrupted;
  const double remaining = time_limit - secondsSince(start);
  if (remaining <= 0.0) return Status::kTimeLimit;

  detail::ShotSolveFunction solve = shotBackend(last_error_);
  if (!solve) return Status::kError;

  detail::ShotOptions options;
  options.model_name = model_name_;
  options.output_flag = output_flag_;
  options.threads = threads_;
  options.mip_rel_gap = mip_rel_gap_;
  options.mip_abs_gap = mip_abs_gap_;
  options.time_limit = remaining;
  options.terminate_flag = terminate_flag_;

  detail::ShotResult result;
  solve(a, options, result);
  if (!result.error.empty()) last_error_ = result.error;
  if (result.status == Status::kOptimal) {
    storeSolution(a, result.col_values);
    obj_val_ = objectiveValue(a, result.col_values);
  }
  return result.status;
}

}  // namespace hs
