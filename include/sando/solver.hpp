/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

// Interface of the trajectory optimization backends, selected at runtime with the "solver"
// parameter: SolverHighsShot (HiGHS + SHOT, "highs_shot") and SolverGurobi ("gurobi"). Each
// backend is compiled when its solver is found at build time (SANDO_HAS_HIGHS_SHOT,
// SANDO_HAS_GUROBI). Callers only use this header, so they do not depend on either solver's headers.

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <decomp_rviz_plugins/data_ros_utils.hpp>
#include <sando/sando_type.hpp>
#include "hgp/termcolor.hpp"
#include "timer.hpp"

using namespace termcolor;
typedef timer::Timer MyTimer;

class TrajSolver {
 public:
  // Sub-step timing breakdown from last generateNewTrajectory call (ms)
  struct SolveTimingBreakdown {
    double findDT_ms{0.0};
    double setX_ms{0.0};
    double polytopes_ms{0.0};
    double dynamic_ms{0.0};
    double objective_ms{0.0};
    double mapsize_ms{0.0};
    double callOptimizer_ms{0.0};
    double postsolve_ms{0.0};
  };
  SolveTimingBreakdown last_solve_timing_;

  virtual ~TrajSolver() = default;

  /** @brief Set the planner name to select the formulation (e.g., "SANDO" or "FASTER"). */
  virtual void setPlannerName(const std::string& name) = 0;

  /** @brief Set the number of threads the solver uses for one optimization. */
  virtual void setGurobiThreads(int num_threads) = 0;

  /** @brief Initialize the solver with planner parameters (segment count, limits, bounds, etc.). */
  virtual void initializeSolver(const Parameters& par) = 0;

  /** @brief Request the running optimization to stop (thread-safe). */
  virtual void stopExecution() = 0;

  /** @brief Clear a stop request so the next optimization runs normally. */
  virtual void resetToNominalState() = 0;

  virtual void setT0(double t0) = 0;
  virtual void setX0(const RobotState& data) = 0;
  virtual void setXf(const RobotState& data) = 0;
  virtual double getInitialDt() = 0;
  virtual void setInitialDt(double initial_dt) = 0;
  virtual void setPolytopes(std::vector<LinearConstraint3D> polytopes) = 0;
  virtual void setPolytopesTimeLayered(
      const std::vector<std::vector<LinearConstraint3D>>& polytopes_by_time) = 0;

  /** @brief Build and solve the trajectory problem with time allocation `factor`.
   *  @param solver_error_detected Set when the solver failed (the solver should then be rebuilt).
   *  @param solver_computation_time Output solve time [ms].
   *  @return True if a feasible trajectory was found.
   */
  virtual bool generateNewTrajectory(
      bool& solver_error_detected,
      double& solver_computation_time,
      double factor,
      bool use_single_thread = false) = 0;

  virtual double getFactorThatWorked() = 0;
  virtual double getObjectiveValue() const = 0;
  virtual void getTotalTrajTime(double& total_traj_time) = 0;
  virtual void fillGoalSetPoints() = 0;
  virtual void getGoalSetpoints(std::vector<RobotState>& goal_setpoints) = 0;
  virtual void getPieceWisePol(PieceWisePol& pwp) = 0;
  virtual void getControlPoints(std::vector<Eigen::Matrix<double, 3, 4>>& cps) = 0;
};

// Kept for the callers written against the build-time selection
using SandoSolver = TrajSolver;

/** @brief Create the backend named `name`: "gurobi", "highs_shot", or "auto" (or empty) for the
 *  first one built, in that order. Throws std::runtime_error if it is unknown, not built, or cannot
 *  be created (e.g., no Gurobi license).
 */
std::shared_ptr<TrajSolver> createTrajSolver(const std::string& name);

/** @brief Resolve "auto"/empty to the backend createTrajSolver() would use; other names unchanged. */
std::string resolveTrajSolverName(const std::string& name);

/** @brief Names of the backends compiled into this build, in the order "auto" prefers them,
 *  e.g. {"gurobi", "highs_shot"}.
 */
std::vector<std::string> availableTrajSolvers();
