/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

// Gazebo Sim (gz) version of dynamic_obstacles_world_plugin.cpp: a single world system spawns
// every dynamic obstacle listed in a JSON file and moves them along their traj_x/y/z expressions.
//
// Usage (inside <world>):
//   <plugin filename="sando_gz_dynamic_obstacles" name="sando::DynamicObstacles">
//     <json_path>/tmp/sando_obstacles.json</json_path>
//     <time_source>sim</time_source>  <!-- "sim" (default) or "wall" -->
//   </plugin>

#include <chrono>
#include <fstream>
#include <list>
#include <memory>
#include <sstream>
#include <string>

#include <gz/common/Console.hh>
#include <gz/math/Pose3.hh>
#include <gz/plugin/Register.hh>
#include <gz/sim/EntityComponentManager.hh>
#include <gz/sim/EventManager.hh>
#include <gz/sim/SdfEntityCreator.hh>
#include <gz/sim/System.hh>
#include <gz/sim/components/Pose.hh>
#include <nlohmann/json.hpp>
#include <sdf/Root.hh>
#include <sim/exprtk.hpp>

using json = nlohmann::json;

namespace sando {

// ------------------------------------------------------------------ obstacle
struct DynObs {
  std::string name;
  double x0{0}, y0{0}, z0{0};
  double bbox_x{0.8}, bbox_y{0.8}, bbox_z{0.8};

  std::string traj_x_str, traj_y_str, traj_z_str;

  // exprtk
  double t_var{0.0};
  exprtk::symbol_table<double> sym;
  exprtk::expression<double> ex, ey, ez;
  bool compiled{false};

  // gz entity (set after spawn)
  gz::sim::Entity entity{gz::sim::kNullEntity};

  void compile() {
    sym.clear();
    sym.add_variable("t", t_var);
    sym.add_constants();

    ex.register_symbol_table(sym);
    ey.register_symbol_table(sym);
    ez.register_symbol_table(sym);

    exprtk::parser<double> p;
    if (!p.compile(traj_x_str, ex))
      gzerr << "[DynObsSystem] compile traj_x failed: " << traj_x_str << "\n";
    if (!p.compile(traj_y_str, ey))
      gzerr << "[DynObsSystem] compile traj_y failed: " << traj_y_str << "\n";
    if (!p.compile(traj_z_str, ez))
      gzerr << "[DynObsSystem] compile traj_z failed: " << traj_z_str << "\n";

    compiled = true;
  }

  inline void evaluate(double t, double& x, double& y, double& z) {
    if (!compiled) {
      x = x0;
      y = y0;
      z = z0;
      return;
    }
    t_var = t;
    x = ex.value();
    y = ey.value();
    z = ez.value();
  }
};

// ------------------------------------------------------------------ system
class DynamicObstacles : public gz::sim::System,
                         public gz::sim::ISystemConfigure,
                         public gz::sim::ISystemPreUpdate {
 public:
  void Configure(const gz::sim::Entity& entity, const std::shared_ptr<const sdf::Element>& sdf,
                 gz::sim::EntityComponentManager& /*ecm*/,
                 gz::sim::EventManager& event_mgr) override {
    world_ = entity;
    event_mgr_ = &event_mgr;

    if (sdf->HasElement("json_path")) json_path_ = sdf->Get<std::string>("json_path");
    if (sdf->HasElement("time_source")) use_wall_time_ = sdf->Get<std::string>("time_source") == "wall";

    gzmsg << "[DynObsSystem] Will read obstacles from: " << json_path_ << " (time source: "
          << (use_wall_time_ ? "wall" : "sim") << ")\n";
  }

  void PreUpdate(const gz::sim::UpdateInfo& info, gz::sim::EntityComponentManager& ecm) override {
    // Phase 1: wait for the JSON file (it may be written after gz starts)
    if (!json_loaded_) {
      if (poll_counter_++ % 100 == 0 && LoadJSON()) SpawnModels(ecm);
      return;
    }

    // Phase 2: move every obstacle.
    // Sim time is the clock that drives the sensors; dynamic_forest_node must then run with
    // use_sim_time:=true (and /clock bridged). Wall time matches nodes running on the system clock.
    double t;
    if (use_wall_time_) {
      t = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
    } else {
      t = std::chrono::duration<double>(info.simTime).count();
    }

    for (auto& o : obstacles_) {
      double x, y, z;
      o.evaluate(t, x, y, z);
      // The models are static with no collisions, so the pose is written directly instead of
      // going through the physics system
      auto* pose = ecm.Component<gz::sim::components::Pose>(o.entity);
      if (!pose) continue;
      pose->Data() = gz::math::Pose3d(x, y, z, 0, 0, 0);
      ecm.SetChanged(o.entity, gz::sim::components::Pose::typeId,
                     gz::sim::ComponentState::PeriodicChange);
    }
  }

 private:
  // ---- helpers ----
  bool LoadJSON() {
    std::ifstream ifs(json_path_);
    if (!ifs.is_open()) {
      gzmsg << "[DynObsSystem] Cannot open: " << json_path_ << "\n";
      return false;
    }

    try {
      json J = json::parse(ifs);
      // Accept both the shared format {"metadata": ..., "obstacles": [...]} and a flat list
      if (J.is_object() && J.contains("obstacles")) J = J["obstacles"];

      for (auto& item : J) {
        // std::list: exprtk expressions hold raw pointers to t_var, so the elements must not move
        DynObs& o = obstacles_.emplace_back();
        o.name = item.value("name", "");
        o.x0 = item.value("x0", 0.0);
        o.y0 = item.value("y0", 0.0);
        o.z0 = item.value("z0", 0.0);
        o.bbox_x = item.value("size_x", item.value("size", 0.8));
        o.bbox_y = item.value("size_y", item.value("size", 0.8));
        o.bbox_z = item.value("size_z", item.value("size", 0.8));
        o.traj_x_str = item.value("traj_x", std::to_string(o.x0));
        o.traj_y_str = item.value("traj_y", std::to_string(o.y0));
        o.traj_z_str = item.value("traj_z", std::to_string(o.z0));
        o.compile();
      }

      gzmsg << "[DynObsSystem] Loaded " << obstacles_.size() << " obstacles from " << json_path_
            << "\n";
      json_loaded_ = true;
      return true;
    } catch (const std::exception& e) {
      gzerr << "[DynObsSystem] JSON parse error: " << e.what() << "\n";
      obstacles_.clear();
      return false;
    }
  }

  void SpawnModels(gz::sim::EntityComponentManager& ecm) {
    gz::sim::SdfEntityCreator creator(ecm, *event_mgr_);

    for (auto& o : obstacles_) {
      std::ostringstream ss;
      // No <collision> element: obstacles pass through each other and the drone.
      // The D435 depth camera renders <visual> geometry, so the point cloud still sees them.
      ss << "<?xml version='1.0' ?>"
         << "<sdf version='1.9'>"
         << "<model name='" << o.name << "'>"
         << "  <static>true</static>"
         << "  <pose>" << o.x0 << " " << o.y0 << " " << o.z0 << " 0 0 0</pose>"
         << "  <link name='link'>"
         << "    <visual name='vis'>"
         << "      <geometry><box><size>" << o.bbox_x << " " << o.bbox_y << " " << o.bbox_z
         << "</size></box></geometry>"
         << "      <material>"
         << "        <ambient>0.5 0.5 0.5 1</ambient>"
         << "        <diffuse>0.5 0.5 0.5 1</diffuse>"
         << "      </material>"
         << "    </visual>"
         << "  </link>"
         << "</model>"
         << "</sdf>";

      sdf::Root root;
      auto errors = root.LoadSdfString(ss.str());
      if (!errors.empty() || !root.Model()) {
        gzerr << "[DynObsSystem] Invalid SDF for obstacle " << o.name << "\n";
        continue;
      }
      o.entity = creator.CreateEntities(root.Model());
      creator.SetParent(o.entity, world_);
    }

    gzmsg << "[DynObsSystem] Spawned " << obstacles_.size() << " box models\n";
  }

  // ---- state ----
  gz::sim::Entity world_{gz::sim::kNullEntity};
  gz::sim::EventManager* event_mgr_{nullptr};
  std::list<DynObs> obstacles_;
  std::string json_path_{"/tmp/sando_obstacles.json"};
  bool use_wall_time_{false};
  bool json_loaded_{false};
  int poll_counter_{0};
};

}  // namespace sando

GZ_ADD_PLUGIN(sando::DynamicObstacles, gz::sim::System, sando::DynamicObstacles::ISystemConfigure,
              sando::DynamicObstacles::ISystemPreUpdate)
GZ_ADD_PLUGIN_ALIAS(sando::DynamicObstacles, "sando::DynamicObstacles")
