/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

// Selects the trajectory optimization backend at build time. SolverHighsShot (HiGHS + SHOT) is the
// default; configuring with -DSANDO_USE_GUROBI=ON defines SANDO_USE_GUROBI and selects
// SolverGurobi. Both classes have the same interface.

#pragma once

#ifdef SANDO_USE_GUROBI
#include <sando/gurobi_solver.hpp>
using SandoSolver = SolverGurobi;
#else
#include <sando/highs_shot_solver.hpp>
using SandoSolver = SolverHighsShot;
#endif
