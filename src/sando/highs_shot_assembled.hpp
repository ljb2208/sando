/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

// Internal to the HiGHS/SHOT modeling layer: the model assembled for the solver backends, and the
// interface of the SHOT backend plugin (libsando_shot_backend.so).
//
// SHOT is only used through that plugin, which is loaded at runtime with
// dlopen(RTLD_LOCAL | RTLD_DEEPBIND) instead of being linked. SHOT compiles in its own copy of
// spdlog and exports its symbols; linked into a ROS node, ROS' logging (which uses the system
// spdlog, a different version) binds to SHOT's copy and the node aborts in rclcpp::init. Loaded
// this way, SHOT's symbols stay private to the plugin and SHOT binds to its own spdlog.

#pragma once

#include <atomic>
#include <string>
#include <vector>
#include <sando/highs_shot_model.hpp>

namespace hs {
namespace detail {

// A sparse affine expression over model columns: sum vals[k] * x[cols[k]] + constant.
struct SparseAffine {
  std::vector<int> cols;
  std::vector<double> vals;
  double constant{0.0};
};

// The term coeff * expr^2 of a sum of squares.
struct SquareTerm {
  double coeff;
  SparseAffine expr;
};

// Row bounds [lo, up] for "activity + constant sense 0".
inline void senseToBounds(char sense, double constant, double& lo, double& up) {
  lo = -kInf;
  up = kInf;
  if (sense == kLessEqual || sense == kEqual) up = -constant;
  if (sense == kGreaterEqual || sense == kEqual) lo = -constant;
}

// A model assembled for the backends, in terms of the active (column-indexed) variables.
struct AssembledModel {
  int num_cols{0};
  std::vector<int> col_slot;  // column -> variable slot
  std::vector<int> slot_col;  // variable slot -> column (-1 if inactive)
  std::vector<double> lb, ub, cost;
  std::vector<char> is_binary;
  bool has_integer{false};  // some binary variable is not fixed
  int sense{kMinimize};
  double offset{0.0};

  // Linear rows in compressed row format
  std::vector<double> row_lo, row_up;
  std::vector<int> row_start{0};
  std::vector<int> row_index;
  std::vector<double> row_value;

  // Quadratic objective terms coeff * x_i * x_j with i <= j, merged
  std::vector<int> obj_qi, obj_qj;
  std::vector<double> obj_qv;

  // The objective as cost' x + offset + sum of squares, when it has that form (used for SHOT)
  bool obj_is_sos{false};
  std::vector<double> obj_sos_cost;
  double obj_sos_offset{0.0};
  std::vector<SquareTerm> obj_squares;

  // Quadratic constraints "lin + quad sense 0", expanded (and as linear + sum of squares if possible)
  struct QRow {
    char sense;
    SparseAffine lin;
    std::vector<int> qi, qj;
    std::vector<double> qv;
    bool is_sos{false};
    SparseAffine sos_lin;
    std::vector<SquareTerm> squares;
  };
  std::vector<QRow> qrows;

  int numRows() const { return static_cast<int>(row_lo.size()); }

  void addRow(const std::vector<int>& cols, const std::vector<double>& vals, double lo, double up) {
    row_index.insert(row_index.end(), cols.begin(), cols.end());
    row_value.insert(row_value.end(), vals.begin(), vals.end());
    row_start.push_back(static_cast<int>(row_index.size()));
    row_lo.push_back(lo);
    row_up.push_back(up);
  }
};

// ---- SHOT backend plugin interface ----

struct ShotOptions {
  std::string model_name;
  bool output_flag{false};
  int threads{0};
  double mip_rel_gap{1e-4};
  double mip_abs_gap{1e-6};
  double time_limit{kInf};  // seconds
  const std::atomic<bool>* terminate_flag{nullptr};
};

struct ShotResult {
  Status status{Status::kError};
  std::vector<double> col_values;  // indexed by column, set when status is kOptimal
  std::string error;
};

// Signature of the plugin's entry point, exported as kShotSolveSymbol. The version in the symbol
// name changes whenever AssembledModel, ShotOptions or ShotResult change.
using ShotSolveFunction = void (*)(const AssembledModel&, const ShotOptions&, ShotResult&);
constexpr const char* kShotSolveSymbol = "sando_hs_shot_solve_v1";

}  // namespace detail
}  // namespace hs
