/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

// Hybrid A* for Ackermann-steered vehicles: searches over (x, y, heading) with bicycle-model arcs
// (steering in [-max_steer, max_steer], forward and optionally reverse), so the path respects the
// minimum turning radius. The vehicle's rectangular footprint is checked against the map at every
// sub-step. The heuristic is an obstacle-aware 2D distance (Dijkstra from the goal), which keeps
// the expansions low in clutter. Independent of ROS: the map is given as query functions.

#pragma once

#include <functional>
#include <string>
#include <vector>

namespace sando {

struct HybridAStarConfig {
  double wheelbase = 0.5;         // [m]
  double max_steer = 0.52;        // [rad] max front-wheel angle
  double step = 0.5;              // [m] arc length of one motion primitive
  int heading_bins = 72;          // heading resolution of the closed set (5 deg)
  bool allow_reverse = false;     // reverse primitives
  double reverse_penalty = 2.0;   // cost factor for distance driven backward
  double switch_penalty = 3.0;    // [m] cost of changing between forward and reverse
  double steer_penalty = 0.3;     // cost factor for |steer| / max_steer per metre
  double steer_change_penalty = 0.2;  // cost factor for steering changes per metre
  double unknown_penalty = 1.0;   // extra cost factor for unknown cells
  double heuristic_weight = 1.2;  // > 1 trades optimality for speed
  int max_expansions = 8000;
  double goal_tolerance = 0.6;    // [m]
  double length = 0.8;            // [m] footprint, centred on the vehicle's reference point
  double width = 0.5;             // [m]
};

/** @brief One pose of the path; `reverse` is the direction used to reach it. */
struct HybridAStarPose {
  double x;
  double y;
  double yaw;  // vehicle heading [rad]
  bool reverse;
};

/** @brief Planar map: a grid of nx * ny cells of size res from (x_min, y_min); cells outside it
 *  are blocked. */
struct HybridAStarMap {
  double x_min = 0.0;
  double y_min = 0.0;
  double res = 0.3;
  int nx = 0;
  int ny = 0;
  std::function<bool(double, double)> occupied;  // obstacle at a point (with safety margin)
  std::function<bool(double, double)> unknown;   // unobserved (traversable at extra cost)
};

class HybridAStar {
 public:
  explicit HybridAStar(const HybridAStarConfig& config);

  /** @brief Plan from (sx, sy, syaw) toward (gx, gy) (any final heading).
   *  @param start_dir Direction of the first move: +1 forward, -1 reverse, 0 either (a vehicle
   *                   already moving cannot switch direction without stopping).
   *  @param path Output poses, from the start to the end.
   *  @param reached_goal True if the path ends within goal_tolerance of the goal; false for a
   *                      partial path toward it (when the goal cannot be reached in the budget).
   *  @param info Short description of the result (for logs).
   *  @return True if a path (complete or partial, making progress) was found.
   */
  bool plan(
      const HybridAStarMap& map,
      double sx,
      double sy,
      double syaw,
      int start_dir,
      double gx,
      double gy,
      std::vector<HybridAStarPose>& path,
      bool& reached_goal,
      std::string& info);

  /** @brief Whether the footprint at (x, y, yaw) is collision-free and inside the map. */
  bool footprintFree(const HybridAStarMap& map, double x, double y, double yaw) const;

 private:
  struct Node {
    double x, y, yaw;
    double g, h;
    int parent;
    bool reverse;
    int steer;  // index into steers_
    int depth;
  };

  void ensureFootprint(double res) const;
  bool computeHeuristic(const HybridAStarMap& map, double& gx, double& gy);
  double heuristicAt(const HybridAStarMap& map, double x, double y) const;
  int cellIndex(const HybridAStarMap& map, double x, double y) const;

  HybridAStarConfig cfg_;
  std::vector<double> steers_;
  mutable std::vector<std::pair<double, double>> footprint_;  // body-frame sample points
  mutable double footprint_res_ = -1.0;                       // map resolution they were made for
  std::vector<float> dist_;                           // 2D distance-to-goal per cell
};

}  // namespace sando
