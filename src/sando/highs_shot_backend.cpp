/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

// SHOT backend plugin (libsando_shot_backend.so) of the HiGHS/SHOT modeling layer. This is the only
// code that includes SHOT headers; it is loaded at runtime by highs_shot_model.cpp, see
// highs_shot_assembled.hpp for why. Built with hidden visibility: only the entry point is exported.

#include "highs_shot_assembled.hpp"

#include <algorithm>
#include <cmath>
#include <memory>

// SHOT
#include "Callback.h"
#include "Enums.h"
#include "Environment.h"
#include "Model/Constraints.h"
#include "Model/ObjectiveFunction.h"
#include "Model/Problem.h"
#include "Model/Terms.h"
#include "Model/Variables.h"
#include "Solver.h"
#include "Structs.h"

// Infinite bounds are represented by IEEE infinity, which -ffast-math assumes away
#ifdef __FAST_MATH__
#error "highs_shot_backend.cpp must be compiled without -ffast-math (add -fno-fast-math)"
#endif

namespace hs {
namespace detail {
namespace {

Status solveWithShot(const AssembledModel& a, const ShotOptions& options, ShotResult& result) {
  auto shouldTerminate = [&options] {
    return options.terminate_flag && options.terminate_flag->load();
  };

  // Solver::hasMIPSolver() has no case for HiGHS, so look it up in the supported list instead
  const auto mip_solvers = SHOT::Solver::getSupportedMIPSolvers();
  if (std::find(mip_solvers.begin(), mip_solvers.end(), SHOT::ES_MIPSolver::Highs) ==
      mip_solvers.end()) {
    result.error = "SHOT was not built with HiGHS support";
    return Status::kError;
  }

  auto solver = std::make_unique<SHOT::Solver>();
  const auto log_level = options.output_flag ? SHOT::E_LogLevel::Info : SHOT::E_LogLevel::Off;
  solver->updateSetting("Output.Console.LogLevel", static_cast<int>(log_level));
  solver->updateSetting("Output.File.LogLevel", static_cast<int>(SHOT::E_LogLevel::Off));
  solver->updateLogLevels();
  solver->updateSetting("Dual.MIP.Solver", static_cast<int>(SHOT::ES_MIPSolver::Highs));
  solver->updateSetting("Dual.MIP.NumberOfThreads", options.threads);
  // The trajectory problems are convex (PSD objective, convex quadratic constraints)
  solver->updateSetting("Model.Convexity.AssumeConvex", true);
  solver->updateSetting("Termination.ObjectiveGap.Relative", options.mip_rel_gap);
  solver->updateSetting("Termination.ObjectiveGap.Absolute", options.mip_abs_gap);
  // The cuts accumulated in the initial LP-relaxation phase can make HiGHS' MIP presolve declare
  // the first MILP infeasible on the trajectory problems, so start directly with MILPs.
  solver->updateSetting("Dual.Relaxation.Use", false);
  // Check linear constraints of primal candidates instead of trusting the MIP solver's values;
  // otherwise integrality tolerance times big-M leaks into corridor constraint violations.
  solver->updateSetting("Primal.Tolerance.TrustLinearConstraintValues", false);
  // No NLP solves (fixed-integer primal heuristic and final polish), so Ipopt is not used: the
  // solution is polished with the binaries fixed afterwards anyway (see Model::polishSolution).
  solver->updateSetting("Primal.FixedInteger.Use", false);
  solver->updateSetting("Primal.PolishSolution.NumberOfPoints", 0);
  if (!std::isinf(options.time_limit))
    solver->updateSetting("Termination.TimeLimit", options.time_limit);

  auto toShot = [](double v) {
    if (v == kInf) return SHOT::SHOT_DBL_MAX;
    if (v == -kInf) return SHOT::SHOT_DBL_MIN;
    return v;
  };

  auto problem = std::make_shared<SHOT::Problem>(solver->getEnvironment());
  problem->name = options.model_name.empty() ? std::string("model") : options.model_name;

  std::vector<SHOT::VariablePtr> variables(a.num_cols);
  for (int j = 0; j < a.num_cols; ++j) {
    const auto type = a.is_binary[j] ? SHOT::E_VariableType::Binary : SHOT::E_VariableType::Real;
    variables[j] = std::make_shared<SHOT::Variable>(
        "x" + std::to_string(j), type, toShot(a.lb[j]), toShot(a.ub[j]));
    problem->add(variables[j]);
  }

  // Sums of squares sum_k c_k (a_k' x + b_k)^2 are passed as sum_k c_k y_k^2 with y_k = a_k' x + b_k.
  // The separable form keeps SHOT's supporting hyperplanes well scaled, whereas the expanded
  // Hessians of the trajectory problems are dense and badly conditioned.
  std::vector<SHOT::VariablePtr> aux_variables;
  std::vector<const SparseAffine*> aux_definitions;
  auto squaresToTerms = [&](const std::vector<SquareTerm>& squares) {
    SHOT::QuadraticTerms terms;
    for (const auto& square : squares) {
      auto y = std::make_shared<SHOT::Variable>(
          "y" + std::to_string(aux_variables.size()), SHOT::E_VariableType::Real, SHOT::SHOT_DBL_MIN,
          SHOT::SHOT_DBL_MAX);
      aux_variables.push_back(y);
      aux_definitions.push_back(&square.expr);
      terms.add(std::make_shared<SHOT::QuadraticTerm>(square.coeff, y, y));
    }
    return terms;
  };
  auto linearTerms = [&](const SparseAffine& expr) {
    SHOT::LinearTerms terms;
    for (size_t k = 0; k < expr.cols.size(); ++k)
      terms.add(std::make_shared<SHOT::LinearTerm>(expr.vals[k], variables[expr.cols[k]]));
    return terms;
  };

  const auto direction = a.sense == kMaximize ? SHOT::E_ObjectiveFunctionDirection::Maximize
                                              : SHOT::E_ObjectiveFunctionDirection::Minimize;
  std::shared_ptr<SHOT::QuadraticObjectiveFunction> objective;
  if (a.obj_is_sos) {
    objective = std::make_shared<SHOT::QuadraticObjectiveFunction>(direction, a.obj_sos_offset);
    for (int j = 0; j < a.num_cols; ++j)
      if (a.obj_sos_cost[j] != 0.0)
        objective->add(std::make_shared<SHOT::LinearTerm>(a.obj_sos_cost[j], variables[j]));
    objective->add(squaresToTerms(a.obj_squares));
  } else {
    objective = std::make_shared<SHOT::QuadraticObjectiveFunction>(direction, a.offset);
    for (int j = 0; j < a.num_cols; ++j)
      if (a.cost[j] != 0.0) objective->add(std::make_shared<SHOT::LinearTerm>(a.cost[j], variables[j]));
    for (size_t k = 0; k < a.obj_qv.size(); ++k)
      objective->add(std::make_shared<SHOT::QuadraticTerm>(
          a.obj_qv[k], variables[a.obj_qi[k]], variables[a.obj_qj[k]]));
  }

  std::vector<std::shared_ptr<SHOT::QuadraticConstraint>> quadratic_constraints;
  for (size_t r = 0; r < a.qrows.size(); ++r) {
    const auto& q = a.qrows[r];
    const SparseAffine& lin = q.is_sos ? q.sos_lin : q.lin;
    double lo, up;
    senseToBounds(q.sense, lin.constant, lo, up);
    SHOT::QuadraticTerms quadratic_terms;
    if (q.is_sos) {
      quadratic_terms = squaresToTerms(q.squares);
    } else {
      for (size_t k = 0; k < q.qv.size(); ++k)
        quadratic_terms.add(
            std::make_shared<SHOT::QuadraticTerm>(q.qv[k], variables[q.qi[k]], variables[q.qj[k]]));
    }
    quadratic_constraints.push_back(std::make_shared<SHOT::QuadraticConstraint>(
        "q" + std::to_string(r), linearTerms(lin), quadratic_terms, toShot(lo), toShot(up)));
  }

  for (auto& y : aux_variables) problem->add(y);
  problem->add(objective);

  for (int r = 0; r < a.numRows(); ++r) {
    SHOT::LinearTerms terms;
    for (int k = a.row_start[r]; k < a.row_start[r + 1]; ++k)
      terms.add(std::make_shared<SHOT::LinearTerm>(a.row_value[k], variables[a.row_index[k]]));
    problem->add(std::make_shared<SHOT::LinearConstraint>(
        "r" + std::to_string(r), terms, toShot(a.row_lo[r]), toShot(a.row_up[r])));
  }

  // y_k - a_k' x = b_k
  for (size_t k = 0; k < aux_variables.size(); ++k) {
    const SparseAffine& def = *aux_definitions[k];
    SHOT::LinearTerms terms;
    terms.add(std::make_shared<SHOT::LinearTerm>(1.0, aux_variables[k]));
    for (size_t i = 0; i < def.cols.size(); ++i)
      terms.add(std::make_shared<SHOT::LinearTerm>(-def.vals[i], variables[def.cols[i]]));
    problem->add(std::make_shared<SHOT::LinearConstraint>(
        "y" + std::to_string(k) + "_def", terms, def.constant, def.constant));
  }

  for (auto& c : quadratic_constraints) problem->add(c);

  problem->finalize();

  if (!solver->setProblem(problem)) {
    result.error = "SHOT could not set up the problem";
    return Status::kError;
  }

  solver->registerCallback<SHOT::TerminationCheckContext>(
      [&shouldTerminate](SHOT::TerminationCheckContext& context) {
        if (shouldTerminate()) context.terminate();
      });

  if (!solver->solveProblem()) {
    result.error = "SHOT could not solve the problem";
    return Status::kError;
  }
  if (shouldTerminate()) return Status::kInterrupted;

  const auto reason = solver->getTerminationReason();
  switch (solver->getModelReturnStatus()) {
    case SHOT::E_ModelReturnStatus::OptimalGlobal: {
      if (!solver->hasPrimalSolution()) return Status::kNumeric;
      const auto solution = solver->getPrimalSolution();
      result.col_values.assign(a.num_cols, 0.0);
      for (int j = 0; j < a.num_cols; ++j) {
        const int index = variables[j]->getIndex();
        if (index < 0 || index >= static_cast<int>(solution.point.size())) {
          result.error = "SHOT returned an incomplete solution";
          return Status::kError;
        }
        result.col_values[j] = solution.point[index];
      }
      return Status::kOptimal;
    }
    case SHOT::E_ModelReturnStatus::InfeasibleGlobal:
    case SHOT::E_ModelReturnStatus::InfeasibleLocal:
      return Status::kInfeasible;
    case SHOT::E_ModelReturnStatus::Unbounded:
    case SHOT::E_ModelReturnStatus::UnboundedNoSolution:
      return Status::kUnbounded;
    case SHOT::E_ModelReturnStatus::ErrorUnknown:
    case SHOT::E_ModelReturnStatus::ErrorNoSolution:
      result.error = "SHOT reported an error";
      return Status::kError;
    default:  // FeasibleSolution or NoSolutionReturned
      if (reason == SHOT::E_TerminationReason::TimeLimit) return Status::kTimeLimit;
      if (reason == SHOT::E_TerminationReason::UserAbort) return Status::kInterrupted;
      return Status::kNumeric;
  }
}

}  // namespace
}  // namespace detail
}  // namespace hs

// Entry point looked up with dlsym (see hs::detail::kShotSolveSymbol)
extern "C" __attribute__((visibility("default"))) void sando_hs_shot_solve_v1(
    const hs::detail::AssembledModel& a, const hs::detail::ShotOptions& options,
    hs::detail::ShotResult& result) {
  result = hs::detail::ShotResult();
  try {
    result.status = hs::detail::solveWithShot(a, options, result);
  } catch (const std::exception& e) {
    result.status = hs::Status::kError;
    result.error = std::string("SHOT: ") + e.what();
  } catch (...) {
    result.status = hs::Status::kError;
    result.error = "SHOT: unknown exception";
  }
  if (result.status != hs::Status::kOptimal) result.col_values.clear();
}
