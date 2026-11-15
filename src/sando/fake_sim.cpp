/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <math.h>
#include <algorithm>
#include <chrono>
#include <mutex>
#include <thread>
#include <Eigen/StdVector>
#include <tf2_ros/buffer.h>
#include <tf2_ros/transform_listener.h>
#include "rclcpp/callback_group.hpp"
#include "rclcpp/executors/multi_threaded_executor.hpp"
#include "rclcpp/rclcpp.hpp"
#include "dynus_interfaces/msg/goal.hpp"
#include "dynus_interfaces/msg/state.hpp"
#include "gazebo_msgs/msg/entity_state.hpp"
#include "gazebo_msgs/srv/set_entity_state.hpp"
#ifdef SANDO_HAS_ROS_GZ
#include "ros_gz_interfaces/msg/entity.hpp"
#include "ros_gz_interfaces/srv/set_entity_pose.hpp"
#endif
#include "geometry_msgs/msg/transform_stamped.hpp"
#include "geometry_msgs/msg/twist.hpp"
#include "nav_msgs/msg/odometry.hpp"
#include "tf2/LinearMath/Quaternion.h"
#include "tf2_ros/transform_broadcaster.h"
#include "visualization_msgs/msg/marker.hpp"
#include <ament_index_cpp/get_package_share_directory.hpp>

using namespace std::chrono_literals;

class FakeSim : public rclcpp::Node {
 public:
  FakeSim() : Node("fake_sim"), br_(this, rclcpp::QoS(10).reliable().durability_volatile()) {
    RCLCPP_INFO(this->get_logger(), "Initializing FakeSim...");

    // Initialize callback groups
    cb_group_me_1_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
    cb_group_re_1_ = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);

    // Declare parameters
    this->declare_parameter<std::vector<double>>("start_pos", {0.0, 0.0, 4.0});
    this->declare_parameter<double>("start_yaw", -1.57);
    this->declare_parameter<bool>("send_state_to_gazebo", true);
    // "classic" (Gazebo Classic, gazebo_ros_state plugin) or "gz" (Gazebo Sim, via ros_gz_bridge)
    this->declare_parameter<std::string>("gazebo_backend", "classic");
    this->declare_parameter<std::string>("gz_set_pose_service", "/world/default/set_pose");
    this->declare_parameter<double>("default_goal_z", 0.3);
    this->declare_parameter<int>("visual_level", 0);

    // New parameters for odometry publishing
    this->declare_parameter<bool>("publish_odom", false);
    this->declare_parameter<std::string>("odom_topic", "odom");
    this->declare_parameter<std::string>("odom_frame_id", "map");
    // If empty, we will set base_frame_id_ = target_frame_ (ns_/base_link)
    this->declare_parameter<std::string>("base_frame_id", "");

    // Vehicle model: "ideal" (the state follows SANDO's setpoints exactly) or "ackermann"
    // (kinematic bicycle driven by cmd_vel from ackermann_controller)
    this->declare_parameter<std::string>("vehicle_model", "ideal");
    this->declare_parameter<double>("ackermann_wheelbase", 0.5);
    this->declare_parameter<double>("ackermann_max_steering_deg", 30.0);
    this->declare_parameter<double>("ackermann_max_accel", 2.0);
    this->declare_parameter<std::vector<double>>("rover_size", {0.8, 0.5, 0.3});

    // Get parameters
    auto start_pos = this->get_parameter("start_pos").as_double_array();
    double yaw = this->get_parameter("start_yaw").as_double();
    send_state_to_gazebo_ = this->get_parameter("send_state_to_gazebo").as_bool();
    gazebo_backend_ = this->get_parameter("gazebo_backend").as_string();
    std::string gz_set_pose_service = this->get_parameter("gz_set_pose_service").as_string();
    default_goal_z_ = this->get_parameter("default_goal_z").as_double();
    int visual_level = this->get_parameter("visual_level").as_int();

    publish_odom_ = this->get_parameter("publish_odom").as_bool();
    odom_topic_ = this->get_parameter("odom_topic").as_string();
    odom_frame_id_ = this->get_parameter("odom_frame_id").as_string();
    base_frame_id_param_ = this->get_parameter("base_frame_id").as_string();
    ackermann_ = this->get_parameter("vehicle_model").as_string() == "ackermann";
    wheelbase_ = this->get_parameter("ackermann_wheelbase").as_double();
    max_steer_ = this->get_parameter("ackermann_max_steering_deg").as_double() * M_PI / 180.0;
    max_accel_ = this->get_parameter("ackermann_max_accel").as_double();
    rover_size_ = this->get_parameter("rover_size").as_double_array();
    rover_yaw_ = yaw;

    // Print parameters
    RCLCPP_INFO(
        this->get_logger(), "Start position: %f, %f, %f", start_pos[0], start_pos[1], start_pos[2]);
    RCLCPP_INFO(this->get_logger(), "Start yaw: %f", yaw);
    RCLCPP_INFO(this->get_logger(), "Send state to Gazebo: %d", send_state_to_gazebo_);
    RCLCPP_INFO(this->get_logger(), "Gazebo backend: %s", gazebo_backend_.c_str());
    RCLCPP_INFO(this->get_logger(), "Default goal z: %f", default_goal_z_);
    RCLCPP_INFO(this->get_logger(), "Visual level: %d", visual_level);
    RCLCPP_INFO(this->get_logger(), "Publish odom: %d", publish_odom_);
    RCLCPP_INFO(this->get_logger(), "Odom topic: %s", odom_topic_.c_str());
    RCLCPP_INFO(this->get_logger(), "Odom frame: %s", odom_frame_id_.c_str());
    RCLCPP_INFO(this->get_logger(), "Base frame param: %s", base_frame_id_param_.c_str());

    // Initialize state
    state_ = dynus_interfaces::msg::State();
    state_.header.frame_id = "map";
    state_.pos.x = start_pos[0];
    state_.pos.y = start_pos[1];
    state_.pos.z = start_pos[2];
    double pitch = 0.0, roll = 0.0;
    tf2::Quaternion quat;
    quat.setRPY(roll, pitch, yaw);
    state_.quat.x = quat.x();
    state_.quat.y = quat.y();
    state_.quat.z = quat.z();
    state_.quat.w = quat.w();

    // Get namespace
    ns_ = get_namespace();
    if (!ns_.empty() && ns_[0] == '/') ns_ = ns_.substr(1);

    // Initialize the tf2 buffer and listener
    tf2_buffer_ = std::make_shared<tf2_ros::Buffer>(this->get_clock());
    tf2_listener_ = std::make_shared<tf2_ros::TransformListener>(*tf2_buffer_);

    target_frame_ = ns_ + "/base_link";

    // base frame for odom (child_frame_id)
    if (!base_frame_id_param_.empty())
      base_frame_id_ = base_frame_id_param_;
    else
      base_frame_id_ = target_frame_;

    // Initialize TF message
    t_.header.frame_id = "map";
    t_.child_frame_id = target_frame_;

    // Publishers
    pub_state_ = this->create_publisher<dynus_interfaces::msg::State>(
        "state", rclcpp::QoS(10).reliable().durability_volatile());

    pub_marker_drone_ = this->create_publisher<visualization_msgs::msg::Marker>(
        "drone_marker", rclcpp::QoS(10).reliable().durability_volatile());

    // Optional odometry publisher
    if (publish_odom_) {
      pub_odom_ = this->create_publisher<nav_msgs::msg::Odometry>(
          odom_topic_, rclcpp::QoS(10).reliable().durability_volatile());
      RCLCPP_INFO(
          this->get_logger(), "Odometry publisher created on topic '%s'", odom_topic_.c_str());
    }

    // Subscribers
    sub_goal_ = this->create_subscription<dynus_interfaces::msg::Goal>(
        "goal", 10, std::bind(&FakeSim::goalCallback, this, std::placeholders::_1));
    if (ackermann_) {
      sub_cmd_vel_ = this->create_subscription<geometry_msgs::msg::Twist>(
          "cmd_vel", 10, [this](const geometry_msgs::msg::Twist& m) {
            std::lock_guard<std::mutex> lk(goal_mtx_);
            cmd_v_ = m.linear.x;
            cmd_w_ = m.angular.z;
          });
      RCLCPP_INFO(this->get_logger(), "Vehicle model: ackermann (wheelbase %.2f m, max steering "
                                      "%.1f deg), driven by cmd_vel", wheelbase_,
                  max_steer_ * 180.0 / M_PI);
    }

    // Timer to simulate TF broadcast
    timer_ = this->create_wall_timer(10ms, std::bind(&FakeSim::pubCallback, this), cb_group_me_1_);

    // Gazebo service client (only when needed — actual sending starts in pubCallback
    // after the initial sleep, giving spawn_entity.py time to create the entity)
    if (send_state_to_gazebo_ && gazebo_backend_ == "gz") {
#ifdef SANDO_HAS_ROS_GZ
      gz_client_ = this->create_client<ros_gz_interfaces::srv::SetEntityPose>(gz_set_pose_service);
      while (!gz_client_->wait_for_service(10s)) {
        RCLCPP_INFO(this->get_logger(), "Gazebo service not available, waiting again...");
      }
#else
      RCLCPP_FATAL(this->get_logger(), "gazebo_backend 'gz' requires building with ros_gz_interfaces");
      throw std::runtime_error("fake_sim was built without ros_gz_interfaces");
#endif
    } else if (send_state_to_gazebo_) {
      gazebo_client_ =
          this->create_client<gazebo_msgs::srv::SetEntityState>("/plug/set_entity_state");
      while (!gazebo_client_->wait_for_service(10s)) {
        RCLCPP_INFO(this->get_logger(), "Gazebo service not available, waiting again...");
      }
    }

    // Flag to publish drone marker
    publish_marker_drone_ = (visual_level > 0);

    // Publish the initial state
    std::this_thread::sleep_for(5s);
    state_.header.stamp = this->get_clock()->now();
    pub_state_->publish(state_);

    // Also publish initial odom if enabled
    if (publish_odom_) {
      publishOdometry();
    }

    // Package path
    package_path_ = ament_index_cpp::get_package_share_directory("sando");

    RCLCPP_INFO(this->get_logger(), "Package path: %s", package_path_.c_str());
    RCLCPP_INFO(this->get_logger(), "FakeSim initialized");
  }

 private:
  std::string package_path_;
  std::string ns_;

  rclcpp::CallbackGroup::SharedPtr cb_group_me_1_;
  rclcpp::CallbackGroup::SharedPtr cb_group_re_1_;

  rclcpp::Publisher<dynus_interfaces::msg::State>::SharedPtr pub_state_;
  rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr pub_marker_drone_;
  rclcpp::Publisher<nav_msgs::msg::Odometry>::SharedPtr pub_odom_;

  rclcpp::Subscription<dynus_interfaces::msg::Goal>::SharedPtr sub_goal_;
  rclcpp::Subscription<geometry_msgs::msg::Twist>::SharedPtr sub_cmd_vel_;

  // Ackermann (kinematic bicycle) model
  bool ackermann_{false};
  double wheelbase_{0.5}, max_steer_{0.52}, max_accel_{2.0};
  std::vector<double> rover_size_{0.8, 0.5, 0.3};
  double cmd_v_{0.0}, cmd_w_{0.0};  // latest cmd_vel (guarded by goal_mtx_)
  double rover_v_{0.0}, rover_yaw_{0.0}, rover_steer_{0.0};
  rclcpp::Time last_step_t_{0, 0, RCL_ROS_TIME};
  rclcpp::Client<gazebo_msgs::srv::SetEntityState>::SharedPtr gazebo_client_;
#ifdef SANDO_HAS_ROS_GZ
  rclcpp::Client<ros_gz_interfaces::srv::SetEntityPose>::SharedPtr gz_client_;
#endif
  std::string gazebo_backend_{"classic"};
  rclcpp::TimerBase::SharedPtr timer_;

  dynus_interfaces::msg::State state_;
  bool publish_marker_drone_{false};
  bool send_state_to_gazebo_{true};

  std::shared_ptr<tf2_ros::Buffer> tf2_buffer_;
  std::shared_ptr<tf2_ros::TransformListener> tf2_listener_;
  std::string target_frame_;

  double default_goal_z_{0.3};
  int drone_marker_id_{1};

  tf2_ros::TransformBroadcaster br_;
  geometry_msgs::msg::TransformStamped t_;

  // Odometry options
  bool publish_odom_{false};
  std::string odom_topic_{"odom"};
  std::string odom_frame_id_{"map"};
  std::string base_frame_id_param_{""};
  std::string base_frame_id_{""};

  // Interpolation state: last goal data for extrapolation between planner updates
  std::mutex goal_mtx_;
  bool goal_received_{false};
  rclcpp::Time goal_stamp_;
  // Position, velocity, acceleration at time of last goal
  Eigen::Vector3d goal_pos_{0, 0, 0};
  Eigen::Vector3d goal_vel_{0, 0, 0};
  Eigen::Vector3d goal_acc_{0, 0, 0};
  // Orientation at time of last goal
  tf2::Quaternion goal_quat_{0, 0, 0, 1};

  // This is for ground robot
  void updateStateFromTF() {
    try {
      geometry_msgs::msg::TransformStamped transform =
          tf2_buffer_->lookupTransform("map", target_frame_, tf2::TimePointZero);

      state_.header.stamp = this->get_clock()->now();
      state_.pos.x = transform.transform.translation.x;
      state_.pos.y = transform.transform.translation.y;
      state_.pos.z = default_goal_z_;

      state_.quat.x = transform.transform.rotation.x;
      state_.quat.y = transform.transform.rotation.y;
      state_.quat.z = transform.transform.rotation.z;
      state_.quat.w = transform.transform.rotation.w;

      pub_state_->publish(state_);
    } catch (const tf2::TransformException& ex) {
      RCLCPP_WARN(this->get_logger(), "Could not transform: %s", ex.what());
    }
  }

  void goalCallback(const dynus_interfaces::msg::Goal::SharedPtr data) {
    if (ackermann_) return;  // driven by cmd_vel instead

    // Hopf fibration approach
    Eigen::Vector3d thrust;
    thrust << data->a.x, data->a.y, data->a.z + 9.81;
    Eigen::Vector3d thrust_normaized = thrust.normalized();

    double a = thrust_normaized.x();
    double b = thrust_normaized.y();
    double c = thrust_normaized.z();

    tf2::Quaternion qabc;
    tf2::Quaternion qpsi;

    double tmp = 1 / std::sqrt(2 * (1 + c));

    qabc.setValue(-b * tmp, a * tmp, 0.0, tmp * (1 + c));

    qpsi.setValue(0.0, 0.0, sin(data->yaw / 2), cos(data->yaw / 2));

    tf2::Quaternion w_q_b = qabc * qpsi;

    // Store goal data for interpolation
    {
      std::lock_guard<std::mutex> lk(goal_mtx_);
      goal_stamp_ = this->get_clock()->now();
      goal_pos_ = Eigen::Vector3d(data->p.x, data->p.y, data->p.z);
      goal_vel_ = Eigen::Vector3d(data->v.x, data->v.y, data->v.z);
      goal_acc_ = Eigen::Vector3d(data->a.x, data->a.y, data->a.z);
      goal_quat_ = w_q_b;
      goal_received_ = true;
    }

    // Also update state_ directly (keeps state_.vel current for odom, etc.)
    state_.header.stamp = this->get_clock()->now();
    state_.pos = data->p;
    state_.vel = data->v;
    state_.quat.w = w_q_b.w();
    state_.quat.x = w_q_b.x();
    state_.quat.y = w_q_b.y();
    state_.quat.z = w_q_b.z();
  }

  void sendGazeboState() {
#ifdef SANDO_HAS_ROS_GZ
    if (gz_client_) {
      sendGzState();
      return;
    }
#endif
    auto request = std::make_shared<gazebo_msgs::srv::SetEntityState::Request>();
    request->state.name = ns_;

    request->state.pose.position.x = state_.pos.x;
    request->state.pose.position.y = state_.pos.y;
    request->state.pose.position.z = state_.pos.z;
    request->state.pose.orientation.x = state_.quat.x;
    request->state.pose.orientation.y = state_.quat.y;
    request->state.pose.orientation.z = state_.quat.z;
    request->state.pose.orientation.w = state_.quat.w;

    auto future = gazebo_client_->async_send_request(request);

    try {
      (void)future.get();
    } catch (const std::exception&) {
      // Keep silent (same behavior as your original code)
    }
  }

#ifdef SANDO_HAS_ROS_GZ
  void sendGzState() {
    auto request = std::make_shared<ros_gz_interfaces::srv::SetEntityPose::Request>();
    request->entity.name = ns_;
    request->entity.type = ros_gz_interfaces::msg::Entity::MODEL;

    request->pose.position.x = state_.pos.x;
    request->pose.position.y = state_.pos.y;
    request->pose.position.z = state_.pos.z;
    request->pose.orientation.x = state_.quat.x;
    request->pose.orientation.y = state_.quat.y;
    request->pose.orientation.z = state_.quat.z;
    request->pose.orientation.w = state_.quat.w;

    auto future = gz_client_->async_send_request(request);

    try {
      (void)future.get();
    } catch (const std::exception&) {
      // Keep silent (same behavior as the Gazebo Classic path)
    }
  }
#endif

  void getInterpolatedTransform() {
    t_.header.stamp = this->get_clock()->now();

    std::lock_guard<std::mutex> lk(goal_mtx_);
    if (!goal_received_) {
      // No goal yet, use current state_ as-is
      t_.transform.translation.x = state_.pos.x;
      t_.transform.translation.y = state_.pos.y;
      t_.transform.translation.z = state_.pos.z;
      t_.transform.rotation.x = state_.quat.x;
      t_.transform.rotation.y = state_.quat.y;
      t_.transform.rotation.z = state_.quat.z;
      t_.transform.rotation.w = state_.quat.w;
      return;
    }

    // Compute elapsed time since last goal, clamped to avoid runaway extrapolation
    double dt = (this->get_clock()->now() - goal_stamp_).seconds();
    dt = std::clamp(dt, 0.0, 0.1);  // cap at 100ms (safe extrapolation window)

    // Extrapolate position: p = p0 + v*dt + 0.5*a*dt^2
    Eigen::Vector3d p_interp = goal_pos_ + goal_vel_ * dt + 0.5 * goal_acc_ * dt * dt;

    t_.transform.translation.x = p_interp.x();
    t_.transform.translation.y = p_interp.y();
    t_.transform.translation.z = p_interp.z();

    t_.transform.rotation.x = goal_quat_.x();
    t_.transform.rotation.y = goal_quat_.y();
    t_.transform.rotation.z = goal_quat_.z();
    t_.transform.rotation.w = goal_quat_.w();
  }

  void publishOdometry() {
    if (!publish_odom_ || !pub_odom_) {
      return;
    }

    nav_msgs::msg::Odometry odom;
    odom.header.stamp = this->get_clock()->now();
    odom.header.frame_id = odom_frame_id_;
    odom.child_frame_id = base_frame_id_;

    // Pose from state_
    odom.pose.pose.position.x = state_.pos.x;
    odom.pose.pose.position.y = state_.pos.y;
    odom.pose.pose.position.z = state_.pos.z;

    odom.pose.pose.orientation.x = state_.quat.x;
    odom.pose.pose.orientation.y = state_.quat.y;
    odom.pose.pose.orientation.z = state_.quat.z;
    odom.pose.pose.orientation.w = state_.quat.w;

    // Twist from state_ (if your dynus_interfaces::msg::State vel is in map/world frame,
    // then this is consistent with header.frame_id = odom_frame_id_. If it's body-frame,
    // you may want to rotate it.)
    odom.twist.twist.linear.x = state_.vel.x;
    odom.twist.twist.linear.y = state_.vel.y;
    odom.twist.twist.linear.z = state_.vel.z;

    // Angular velocity is unknown here; leave as zeros.

    pub_odom_->publish(odom);
  }

  /** Integrate the kinematic bicycle model with the latest cmd_vel (speed, yaw rate). */
  void stepAckermann() {
    const rclcpp::Time t = this->get_clock()->now();
    if (last_step_t_.nanoseconds() == 0 || last_step_t_.get_clock_type() != t.get_clock_type()) {
      last_step_t_ = t;
      return;
    }
    const double dt = std::clamp((t - last_step_t_).seconds(), 0.0, 0.1);
    last_step_t_ = t;
    double v_cmd, w_cmd;
    {
      std::lock_guard<std::mutex> lk(goal_mtx_);
      v_cmd = cmd_v_;
      w_cmd = cmd_w_;
    }
    // Speed follows the command with limited acceleration; steering comes from the commanded yaw
    // rate (omega = v tan(delta) / L), limited to the steering range
    rover_v_ += std::clamp(v_cmd - rover_v_, -max_accel_ * dt, max_accel_ * dt);
    if (std::fabs(v_cmd) > 1e-3) rover_steer_ = std::atan(wheelbase_ * w_cmd / v_cmd);
    rover_steer_ = std::clamp(rover_steer_, -max_steer_, max_steer_);
    const double yaw_rate = rover_v_ * std::tan(rover_steer_) / wheelbase_;
    state_.pos.x += rover_v_ * std::cos(rover_yaw_) * dt;
    state_.pos.y += rover_v_ * std::sin(rover_yaw_) * dt;
    rover_yaw_ += yaw_rate * dt;
    state_.vel.x = rover_v_ * std::cos(rover_yaw_);
    state_.vel.y = rover_v_ * std::sin(rover_yaw_);
    state_.vel.z = 0.0;
    tf2::Quaternion q;
    q.setRPY(0.0, 0.0, rover_yaw_);
    state_.quat.x = q.x();
    state_.quat.y = q.y();
    state_.quat.z = q.z();
    state_.quat.w = q.w();
  }

  void pubCallback() {
    if (ackermann_) stepAckermann();

    // Compute interpolated transform (extrapolates between planner updates)
    getInterpolatedTransform();
    br_.sendTransform(t_);

    // Publish odometry (optional)
    if (publish_odom_) {
      publishOdometry();
    }

    // Publish drone marker using interpolated position
    if (publish_marker_drone_) {
      pub_marker_drone_->publish(getDroneMarker());
    }

    // Send the state to Gazebo
    if (send_state_to_gazebo_) {
      std::thread(&FakeSim::sendGazeboState, this).detach();
    }

    // Publish the state
    state_.header.stamp = this->get_clock()->now();
    pub_state_->publish(state_);
  }

  visualization_msgs::msg::Marker getDroneMarker() {
    visualization_msgs::msg::Marker marker;
    marker.id = drone_marker_id_;
    marker.ns = std::string("mesh_") + this->get_namespace();
    marker.header.frame_id = "map";
    marker.header.stamp = this->get_clock()->now();
    marker.type = marker.MESH_RESOURCE;
    marker.action = marker.ADD;

    // Use the interpolated position from the TF transform
    marker.pose.position.x = t_.transform.translation.x;
    marker.pose.position.y = t_.transform.translation.y;
    marker.pose.position.z = t_.transform.translation.z;
    marker.pose.orientation.x = t_.transform.rotation.x;
    marker.pose.orientation.y = t_.transform.rotation.y;
    marker.pose.orientation.z = t_.transform.rotation.z;
    marker.pose.orientation.w = t_.transform.rotation.w;

    if (ackermann_) {  // rover body
      marker.type = marker.CUBE;
      marker.scale.x = rover_size_.size() > 0 ? rover_size_[0] : 0.8;
      marker.scale.y = rover_size_.size() > 1 ? rover_size_[1] : 0.5;
      marker.scale.z = rover_size_.size() > 2 ? rover_size_[2] : 0.3;
      marker.color.r = 0.16f;
      marker.color.g = 0.47f;
      marker.color.b = 0.84f;
      marker.color.a = 1.0f;
      return marker;
    }
    marker.mesh_use_embedded_materials = true;
    marker.mesh_resource = "package://sando/meshes/quadrotor/quadrotor.dae";
    marker.scale.x = 0.75;
    marker.scale.y = 0.75;
    marker.scale.z = 0.75;

    return marker;
  }
};

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);
  auto node = std::make_shared<FakeSim>();
  rclcpp::executors::MultiThreadedExecutor executor;
  executor.add_node(node);
  executor.spin();
  rclcpp::shutdown();
  return 0;
}
