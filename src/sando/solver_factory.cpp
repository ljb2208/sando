/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <sando/solver.hpp>

#include <stdexcept>

#ifdef SANDO_HAS_HIGHS_SHOT
#include <sando/highs_shot_solver.hpp>
#endif
#ifdef SANDO_HAS_GUROBI
#include <sando/gurobi_solver.hpp>
#endif

// In order of preference for "auto": Gurobi solves the trajectory problems much faster
std::vector<std::string> availableTrajSolvers() {
  std::vector<std::string> names;
#ifdef SANDO_HAS_GUROBI
  names.emplace_back("gurobi");
#endif
#ifdef SANDO_HAS_HIGHS_SHOT
  names.emplace_back("highs_shot");
#endif
  return names;
}

std::string resolveTrajSolverName(const std::string& name) {
  if (!name.empty() && name != "auto") return name;
  const auto names = availableTrajSolvers();
  return names.empty() ? std::string() : names.front();
}

std::shared_ptr<TrajSolver> createTrajSolver(const std::string& name) {
  const std::string resolved = resolveTrajSolverName(name);

  std::string available;
  for (const auto& n : availableTrajSolvers()) available += (available.empty() ? "" : ", ") + n;
  if (available.empty()) available = "none";

  if (resolved == "highs_shot" || resolved == "gurobi") {
#ifdef SANDO_HAS_HIGHS_SHOT
    if (resolved == "highs_shot") return std::make_shared<SolverHighsShot>();
#endif
#ifdef SANDO_HAS_GUROBI
    if (resolved == "gurobi") {
      // The model is created with the shared Gurobi environment in the constructor, which throws
      // without a usable license
      try {
        return std::make_shared<SolverGurobi>();
      } catch (const GRBException& e) {
        throw std::runtime_error(
            "could not create the Gurobi solver (" + e.getMessage() +
            "); check the Gurobi license, or use solver:=highs_shot");
      }
    }
#endif
    throw std::runtime_error(
        "solver '" + resolved + "' was not built into this SANDO (available: " + available + ")");
  }
  throw std::runtime_error(
      "unknown solver '" + name + "' (use auto, highs_shot or gurobi; built: " + available + ")");
}
