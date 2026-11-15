/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include "hgp/hybrid_astar.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <sstream>

namespace sando {

namespace {
double wrapAngle(double a) { return std::atan2(std::sin(a), std::cos(a)); }
constexpr float kInf = std::numeric_limits<float>::infinity();
}  // namespace

HybridAStar::HybridAStar(const HybridAStarConfig& config) : cfg_(config) {
  for (double f : {-1.0, -0.5, 0.0, 0.5, 1.0}) steers_.push_back(f * cfg_.max_steer);
}

int HybridAStar::cellIndex(const HybridAStarMap& map, double x, double y) const {
  const int ix = static_cast<int>(std::floor((x - map.x_min) / map.res));
  const int iy = static_cast<int>(std::floor((y - map.y_min) / map.res));
  if (ix < 0 || iy < 0 || ix >= map.nx || iy >= map.ny) return -1;
  return iy * map.nx + ix;
}

void HybridAStar::ensureFootprint(double res) const {
  if (footprint_res_ == res && !footprint_.empty()) return;
  // A grid over the rectangle with samples at most one cell apart
  footprint_.clear();
  const int nu = std::max(2, static_cast<int>(std::ceil(cfg_.length / res)) + 1);
  const int nv = std::max(2, static_cast<int>(std::ceil(cfg_.width / res)) + 1);
  for (int i = 0; i < nu; ++i)
    for (int j = 0; j < nv; ++j)
      footprint_.emplace_back(
          -cfg_.length / 2 + cfg_.length * i / (nu - 1), -cfg_.width / 2 + cfg_.width * j / (nv - 1));
  footprint_res_ = res;
}

bool HybridAStar::footprintFree(const HybridAStarMap& map, double x, double y, double yaw) const {
  ensureFootprint(map.res);
  const double c = std::cos(yaw), s = std::sin(yaw);
  for (const auto& [u, v] : footprint_) {
    const double px = x + u * c - v * s, py = y + u * s + v * c;
    if (cellIndex(map, px, py) < 0 || map.occupied(px, py)) return false;
  }
  return true;
}

// Obstacle-aware 2D distance to the goal (8-connected Dijkstra over the grid). An occupied goal
// moves to the nearest free cell within 2 m.
bool HybridAStar::computeHeuristic(const HybridAStarMap& map, double& gx, double& gy) {
  const int n = map.nx * map.ny;
  std::vector<char> blocked(n), unknown(n);
  for (int iy = 0; iy < map.ny; ++iy)
    for (int ix = 0; ix < map.nx; ++ix) {
      const double x = map.x_min + (ix + 0.5) * map.res, y = map.y_min + (iy + 0.5) * map.res;
      blocked[iy * map.nx + ix] = map.occupied(x, y);
      unknown[iy * map.nx + ix] = map.unknown ? map.unknown(x, y) : 0;
    }

  int goal = cellIndex(map, gx, gy);
  if (goal < 0) return false;
  if (blocked[goal]) {
    const int r = static_cast<int>(std::ceil(2.0 / map.res));
    const int gix = goal % map.nx, giy = goal / map.nx;
    int best = -1;
    double best_d = 1e9;
    for (int dy = -r; dy <= r; ++dy)
      for (int dx = -r; dx <= r; ++dx) {
        const int ix = gix + dx, iy = giy + dy;
        if (ix < 0 || iy < 0 || ix >= map.nx || iy >= map.ny) continue;
        const int k = iy * map.nx + ix;
        const double d = std::hypot(dx, dy);
        if (!blocked[k] && d < best_d) best_d = d, best = k;
      }
    if (best < 0) return false;
    goal = best;
    gx = map.x_min + (goal % map.nx + 0.5) * map.res;
    gy = map.y_min + (goal / map.nx + 0.5) * map.res;
  }

  dist_.assign(n, kInf);
  using QE = std::pair<float, int>;
  std::priority_queue<QE, std::vector<QE>, std::greater<QE>> q;
  dist_[goal] = 0.0f;
  q.push({0.0f, goal});
  const int dxs[8] = {1, -1, 0, 0, 1, 1, -1, -1}, dys[8] = {0, 0, 1, -1, 1, -1, 1, -1};
  while (!q.empty()) {
    auto [d, k] = q.top();
    q.pop();
    if (d > dist_[k]) continue;
    const int ix = k % map.nx, iy = k / map.nx;
    for (int i = 0; i < 8; ++i) {
      const int jx = ix + dxs[i], jy = iy + dys[i];
      if (jx < 0 || jy < 0 || jx >= map.nx || jy >= map.ny) continue;
      const int j = jy * map.nx + jx;
      if (blocked[j]) continue;
      const float step = static_cast<float>(map.res * (i < 4 ? 1.0 : M_SQRT2) *
                                            (unknown[j] ? 1.0 + cfg_.unknown_penalty : 1.0));
      if (d + step < dist_[j]) {
        dist_[j] = d + step;
        q.push({dist_[j], j});
      }
    }
  }
  return true;
}

double HybridAStar::heuristicAt(const HybridAStarMap& map, double x, double y) const {
  const int k = cellIndex(map, x, y);
  return (k < 0) ? kInf : dist_[k];
}

bool HybridAStar::plan(
    const HybridAStarMap& map,
    double sx,
    double sy,
    double syaw,
    int start_dir,
    double gx,
    double gy,
    std::vector<HybridAStarPose>& path,
    bool& reached_goal,
    std::string& info) {
  path.clear();
  reached_goal = false;
  std::ostringstream os;
  if (map.nx <= 0 || map.ny <= 0 || !map.occupied) {
    info = "empty map";
    return false;
  }

  ensureFootprint(map.res);

  if (!computeHeuristic(map, gx, gy)) {
    info = "goal outside the map or enclosed";
    return false;
  }
  const double h_start = heuristicAt(map, sx, sy);
  if (!std::isfinite(h_start)) {
    info = "goal unreachable from the start (2D)";
    return false;
  }

  const int nh = cfg_.heading_bins;
  auto key = [&](double x, double y, double yaw, bool rev) -> int {
    const int c = cellIndex(map, x, y);
    if (c < 0) return -1;
    int h = static_cast<int>(std::floor((wrapAngle(yaw) + M_PI) / (2 * M_PI) * nh)) % nh;
    return (c * nh + h) * 2 + (rev ? 1 : 0);
  };
  std::vector<float> best_g(static_cast<size_t>(map.nx) * map.ny * nh * 2, kInf);
  std::vector<char> closed(best_g.size(), 0);

  std::vector<Node> nodes;
  nodes.reserve(cfg_.max_expansions * 4);
  nodes.push_back({sx, sy, syaw, 0.0, h_start, -1, start_dir < 0, 2, 0});
  using QE = std::pair<double, int>;
  std::priority_queue<QE, std::vector<QE>, std::greater<QE>> open;
  open.push({cfg_.heuristic_weight * h_start, 0});

  const double sub = std::min(map.res * 0.5, cfg_.step);
  const int n_sub = std::max(1, static_cast<int>(std::ceil(cfg_.step / sub)));
  const double ds = cfg_.step / n_sub;
  int best_node = 0, goal_node = -1, expansions = 0;

  while (!open.empty() && expansions < cfg_.max_expansions) {
    const int idx = open.top().second;
    open.pop();
    const Node cur = nodes[idx];
    const int k = key(cur.x, cur.y, cur.yaw, cur.reverse);
    if (k >= 0) {
      if (closed[k]) continue;
      closed[k] = 1;
    }
    ++expansions;
    if (cur.h < nodes[best_node].h) best_node = idx;
    if (std::hypot(cur.x - gx, cur.y - gy) <= cfg_.goal_tolerance) {
      goal_node = idx;
      break;
    }

    for (int dir : {1, -1}) {
      if (dir < 0 && !cfg_.allow_reverse) continue;
      if (cur.depth == 0 && start_dir != 0 && dir != start_dir) continue;
      const bool rev = dir < 0;
      for (size_t si = 0; si < steers_.size(); ++si) {
        const double steer = steers_[si];
        double x = cur.x, y = cur.y, yaw = cur.yaw;
        bool ok = true;
        for (int s = 0; s < n_sub && ok; ++s) {
          // Exact bicycle-model arc of length ds
          const double d = dir * ds;
          if (std::fabs(steer) < 1e-6) {
            x += d * std::cos(yaw);
            y += d * std::sin(yaw);
          } else {
            const double R = cfg_.wheelbase / std::tan(steer);
            const double dyaw = d / R;
            x += R * (std::sin(yaw + dyaw) - std::sin(yaw));
            y += R * (std::cos(yaw) - std::cos(yaw + dyaw));
            yaw = wrapAngle(yaw + dyaw);
          }
          // The first move out of the start only checks the reference point: the inflated map
          // can overlap the footprint where the vehicle already is
          ok = (cur.depth == 0) ? (cellIndex(map, x, y) >= 0 && !map.occupied(x, y))
                                : footprintFree(map, x, y, yaw);
        }
        if (!ok) continue;
        const int nk = key(x, y, yaw, rev);
        if (nk < 0 || closed[nk]) continue;
        const double h = heuristicAt(map, x, y);
        if (!std::isfinite(h)) continue;

        double cost = cfg_.step * (rev ? cfg_.reverse_penalty : 1.0);
        cost += cfg_.step * cfg_.steer_penalty * std::fabs(steer) / cfg_.max_steer;
        if (cur.depth > 0) {
          if (rev != cur.reverse) cost += cfg_.switch_penalty;
          cost += cfg_.step * cfg_.steer_change_penalty * std::fabs(steer - steers_[cur.steer]) /
                  cfg_.max_steer;
        }
        if (map.unknown && map.unknown(x, y)) cost *= 1.0 + cfg_.unknown_penalty;
        const double g = cur.g + cost;
        if (g >= best_g[nk]) continue;
        best_g[nk] = static_cast<float>(g);
        nodes.push_back({x, y, yaw, g, h, idx, rev, static_cast<int>(si), cur.depth + 1});
        open.push({g + cfg_.heuristic_weight * h, static_cast<int>(nodes.size()) - 1});
      }
    }
  }

  int end = goal_node;
  if (end < 0) {
    // Partial path toward the goal, if it gets meaningfully closer
    if (nodes[best_node].h < h_start - 0.5) end = best_node;
  }
  os << "expansions " << expansions << ", nodes " << nodes.size();
  if (end < 0) {
    os << ", no progress toward the goal";
    info = os.str();
    return false;
  }
  for (int i = end; i >= 0; i = nodes[i].parent)
    path.push_back({nodes[i].x, nodes[i].y, nodes[i].yaw, nodes[i].reverse});
  std::reverse(path.begin(), path.end());
  if (path.size() > 1) path.front().reverse = path[1].reverse;
  reached_goal = goal_node >= 0;
  os << (reached_goal ? ", reached the goal" : ", partial path") << " (" << path.size()
     << " poses)";
  info = os.str();
  return path.size() > 1;
}

}  // namespace sando
