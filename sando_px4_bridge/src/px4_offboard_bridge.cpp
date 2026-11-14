/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <sando_px4_bridge/px4_offboard_bridge.hpp>

using px4_msgs::msg::OffboardControlMode;
using px4_msgs::msg::PositionSetpoint;
using px4_msgs::msg::PositionSetpointTriplet;
using px4_msgs::msg::TrajectorySetpoint;
using px4_msgs::msg::VehicleCommand;
using px4_msgs::msg::VehicleLocalPosition;
using px4_msgs::msg::VehicleOdometry;
using px4_msgs::msg::VehicleStatus;

namespace {
constexpr double kEarthRadius = 6371000.0;  // [m], same as PX4's MapProjection
}

Px4OffboardBridge::Px4OffboardBridge()
    : Node("px4_offboard_bridge"),
      pos_ned_(Eigen::Vector3d::Zero()),
      yaw_ned_(0.0),
      nav_state_(0),
      arming_state_(0),
      ref_lat_(0.0),
      ref_lon_(0.0),
      ref_alt_(0.0),
      goal_stamp_(0, 0, RCL_ROS_TIME),
      hold_pos_ned_(Eigen::Vector3d::Zero()),
      hold_yaw_ned_(0.0),
      setpoint_counter_(0),
      odom_initialized_(false),
      ref_initialized_(false),
      goal_initialized_(false),
      engage_sent_(false) {
  // PX4 topic names (check `ros2 topic list`: newer firmware may append a _vN version suffix)
  this->declare_parameter("offboard_control_mode_topic", "/fmu/in/offboard_control_mode");
  this->declare_parameter("trajectory_setpoint_topic", "/fmu/in/trajectory_setpoint");
  this->declare_parameter("vehicle_command_topic", "/fmu/in/vehicle_command");
  this->declare_parameter("vehicle_odometry_topic", "/fmu/out/vehicle_odometry");
  this->declare_parameter("vehicle_local_position_topic", "/fmu/out/vehicle_local_position");
  this->declare_parameter("vehicle_status_topic", "/fmu/out/vehicle_status");
  this->declare_parameter("position_setpoint_triplet_topic", "/fmu/out/position_setpoint_triplet");

  // Behavior
  this->declare_parameter("frame_id", "map");
  this->declare_parameter("publish_state", true);
  this->declare_parameter("auto_engage", false);
  this->declare_parameter("setpoint_rate_hz", 50.0);
  this->declare_parameter("goal_timeout", 0.5);
  this->declare_parameter("term_goal_update_threshold", 0.5);
  this->declare_parameter("target_system", 1);

  frame_id_ = this->get_parameter("frame_id").as_string();
  publish_state_ = this->get_parameter("publish_state").as_bool();
  auto_engage_ = this->get_parameter("auto_engage").as_bool();
  setpoint_rate_hz_ = this->get_parameter("setpoint_rate_hz").as_double();
  goal_timeout_ = this->get_parameter("goal_timeout").as_double();
  term_goal_update_threshold_ = this->get_parameter("term_goal_update_threshold").as_double();
  target_system_ = static_cast<uint8_t>(this->get_parameter("target_system").as_int());

  std::cout << "auto_engage_: " << auto_engage_ << std::endl;
  std::cout << "setpoint_rate_hz_: " << setpoint_rate_hz_ << std::endl;
  std::cout << "goal_timeout_: " << goal_timeout_ << std::endl;

  // PX4 publishes with best-effort reliability
  rmw_qos_profile_t qos_profile = rmw_qos_profile_sensor_data;
  auto px4_qos = rclcpp::QoS(rclcpp::QoSInitialization(qos_profile.history, 5), qos_profile);

  // Publishers and Subscribers (PX4 side)
  pub_offboard_control_mode_ = this->create_publisher<OffboardControlMode>(
      this->get_parameter("offboard_control_mode_topic").as_string(), 10);
  pub_trajectory_setpoint_ = this->create_publisher<TrajectorySetpoint>(
      this->get_parameter("trajectory_setpoint_topic").as_string(), 10);
  pub_vehicle_command_ = this->create_publisher<VehicleCommand>(
      this->get_parameter("vehicle_command_topic").as_string(), 10);
  sub_odometry_ = this->create_subscription<VehicleOdometry>(
      this->get_parameter("vehicle_odometry_topic").as_string(), px4_qos,
      std::bind(&Px4OffboardBridge::odometryCallback, this, std::placeholders::_1));
  sub_local_position_ = this->create_subscription<VehicleLocalPosition>(
      this->get_parameter("vehicle_local_position_topic").as_string(), px4_qos,
      std::bind(&Px4OffboardBridge::localPositionCallback, this, std::placeholders::_1));
  sub_vehicle_status_ = this->create_subscription<VehicleStatus>(
      this->get_parameter("vehicle_status_topic").as_string(), px4_qos,
      std::bind(&Px4OffboardBridge::vehicleStatusCallback, this, std::placeholders::_1));
  sub_setpoint_triplet_ = this->create_subscription<PositionSetpointTriplet>(
      this->get_parameter("position_setpoint_triplet_topic").as_string(), px4_qos,
      std::bind(&Px4OffboardBridge::setpointTripletCallback, this, std::placeholders::_1));

  // Publishers and Subscribers (SANDO side)
  if (publish_state_) pub_state_ = this->create_publisher<dynus_interfaces::msg::State>("state", 10);
  pub_term_goal_ = this->create_publisher<geometry_msgs::msg::PoseStamped>("term_goal", 10);
  sub_goal_ = this->create_subscription<dynus_interfaces::msg::Goal>(
      "goal", 10, std::bind(&Px4OffboardBridge::goalCallback, this, std::placeholders::_1));

  // Timers (PX4 requires the offboard stream at > 2 Hz)
  timer_ = this->create_wall_timer(std::chrono::duration<double>(1.0 / setpoint_rate_hz_),
                                   std::bind(&Px4OffboardBridge::offboardTimerCallback, this));
}

void Px4OffboardBridge::odometryCallback(const VehicleOdometry::SharedPtr msg) {
  if (msg->pose_frame != VehicleOdometry::POSE_FRAME_NED || std::isnan(msg->position[0]) ||
      std::isnan(msg->q[0]))
    return;

  pos_ned_ = Eigen::Vector3d(msg->position[0], msg->position[1], msg->position[2]);
  Eigen::Quaterniond q_ned_frd(msg->q[0], msg->q[1], msg->q[2], msg->q[3]);  // PX4: [w, x, y, z]
  Eigen::Matrix3d R_ned_frd = q_ned_frd.normalized().toRotationMatrix();
  yaw_ned_ = std::atan2(R_ned_frd(1, 0), R_ned_frd(0, 0));
  odom_initialized_ = true;

  if (!publish_state_) return;

  // NED -> ENU (world) and FRD -> FLU (body)
  Eigen::Matrix3d T_enu_ned;
  T_enu_ned << 0, 1, 0, 1, 0, 0, 0, 0, -1;
  Eigen::Matrix3d T_frd_flu = Eigen::Vector3d(1, -1, -1).asDiagonal();
  Eigen::Quaterniond q_enu_flu(T_enu_ned * R_ned_frd * T_frd_flu);

  Eigen::Vector3d vel(msg->velocity[0], msg->velocity[1], msg->velocity[2]);
  if (msg->velocity_frame == VehicleOdometry::VELOCITY_FRAME_BODY_FRD) vel = R_ned_frd * vel;
  if (vel.hasNaN()) vel.setZero();
  Eigen::Vector3d vel_enu = enuNedSwap(vel);
  Eigen::Vector3d pos_enu = enuNedSwap(pos_ned_);

  dynus_interfaces::msg::State state;
  state.header.stamp = this->now();
  state.header.frame_id = frame_id_;
  state.pos.x = pos_enu.x();
  state.pos.y = pos_enu.y();
  state.pos.z = pos_enu.z();
  state.vel.x = vel_enu.x();
  state.vel.y = vel_enu.y();
  state.vel.z = vel_enu.z();
  state.quat.x = q_enu_flu.x();
  state.quat.y = q_enu_flu.y();
  state.quat.z = q_enu_flu.z();
  state.quat.w = q_enu_flu.w();
  pub_state_->publish(state);
}

void Px4OffboardBridge::localPositionCallback(const VehicleLocalPosition::SharedPtr msg) {
  if (!msg->xy_global || !msg->z_global) return;
  ref_lat_ = msg->ref_lat;
  ref_lon_ = msg->ref_lon;
  ref_alt_ = msg->ref_alt;
  ref_initialized_ = true;
}

void Px4OffboardBridge::vehicleStatusCallback(const VehicleStatus::SharedPtr msg) {
  nav_state_ = msg->nav_state;
  arming_state_ = msg->arming_state;
}

void Px4OffboardBridge::setpointTripletCallback(const PositionSetpointTriplet::SharedPtr msg) {
  // The navigator's current setpoint is the active mission waypoint or "Go to" location. Ignore it
  // while in offboard, where the navigator is inactive and only republishes stale setpoints.
  const PositionSetpoint& sp = msg->current;
  if (!ref_initialized_ || !sp.valid || nav_state_ == VehicleStatus::NAVIGATION_STATE_OFFBOARD)
    return;
  if (sp.type != PositionSetpoint::SETPOINT_TYPE_POSITION &&
      sp.type != PositionSetpoint::SETPOINT_TYPE_LOITER)
    return;

  Eigen::Vector3d goal_ned = globalToLocalNed(sp.lat, sp.lon, sp.alt);
  if (goal_ned.hasNaN()) return;
  if (last_term_goal_ned_ &&
      (goal_ned - *last_term_goal_ned_).norm() < term_goal_update_threshold_)
    return;
  last_term_goal_ned_ = goal_ned;

  Eigen::Vector3d goal_enu = enuNedSwap(goal_ned);
  double yaw_enu = std::isfinite(sp.yaw) ? yawEnuNedSwap(sp.yaw) : 0.0;

  geometry_msgs::msg::PoseStamped term_goal;
  term_goal.header.stamp = this->now();
  term_goal.header.frame_id = frame_id_;
  term_goal.pose.position.x = goal_enu.x();
  term_goal.pose.position.y = goal_enu.y();
  term_goal.pose.position.z = goal_enu.z();
  term_goal.pose.orientation.z = std::sin(yaw_enu / 2.0);
  term_goal.pose.orientation.w = std::cos(yaw_enu / 2.0);
  pub_term_goal_->publish(term_goal);

  RCLCPP_INFO(this->get_logger(), "New PX4 waypoint -> term_goal (ENU): [%.2f, %.2f, %.2f]",
              goal_enu.x(), goal_enu.y(), goal_enu.z());
}

void Px4OffboardBridge::goalCallback(const dynus_interfaces::msg::Goal::SharedPtr msg) {
  goal_ = *msg;
  goal_stamp_ = this->now();
  goal_initialized_ = true;
}

void Px4OffboardBridge::offboardTimerCallback() {
  if (!odom_initialized_) return;

  // Until PX4 is in offboard, keep the hover setpoint on the vehicle and drop any stale SANDO goal,
  // so that switching to offboard holds the current position instead of jumping.
  if (nav_state_ != VehicleStatus::NAVIGATION_STATE_OFFBOARD) {
    hold_pos_ned_ = pos_ned_;
    hold_yaw_ned_ = yaw_ned_;
    if (goal_initialized_ && (this->now() - goal_stamp_).seconds() >= goal_timeout_)
      goal_initialized_ = false;
  }

  // Following the PX4 offboard example: stream setpoints for ~1 s, then request offboard and arm
  if (auto_engage_ && !engage_sent_ && goal_initialized_ &&
      setpoint_counter_ >= static_cast<int>(setpoint_rate_hz_)) {
    publishVehicleCommand(VehicleCommand::VEHICLE_CMD_DO_SET_MODE, 1, 6);
    publishVehicleCommand(VehicleCommand::VEHICLE_CMD_COMPONENT_ARM_DISARM, 1.0);
    RCLCPP_INFO(this->get_logger(), "Offboard mode and arm commands sent");
    engage_sent_ = true;
  }

  publishOffboardControlMode();
  publishTrajectorySetpoint();
  setpoint_counter_++;
}

void Px4OffboardBridge::publishOffboardControlMode() {
  OffboardControlMode msg{};
  msg.position = true;
  msg.velocity = true;
  msg.acceleration = true;
  msg.attitude = false;
  msg.body_rate = false;
  msg.timestamp = nowMicros();
  pub_offboard_control_mode_->publish(msg);
}

void Px4OffboardBridge::publishTrajectorySetpoint() {
  const float nan = std::numeric_limits<float>::quiet_NaN();
  TrajectorySetpoint msg{};
  msg.position = {nan, nan, nan};
  msg.velocity = {nan, nan, nan};
  msg.acceleration = {nan, nan, nan};
  msg.jerk = {nan, nan, nan};
  msg.yawspeed = nan;

  bool goal_fresh =
      goal_initialized_ && (this->now() - goal_stamp_).seconds() < goal_timeout_;

  if (goal_initialized_) {
    // Track the SANDO setpoint. Once it goes stale (SANDO stops publishing at the goal on
    // hardware), hold its last position with velocity/acceleration left to PX4.
    Eigen::Vector3d p = enuNedSwap(Eigen::Vector3d(goal_.p.x, goal_.p.y, goal_.p.z));
    msg.position = {static_cast<float>(p.x()), static_cast<float>(p.y()),
                    static_cast<float>(p.z())};
    msg.yaw = static_cast<float>(yawEnuNedSwap(goal_.yaw));

    if (goal_fresh) {
      Eigen::Vector3d v = enuNedSwap(Eigen::Vector3d(goal_.v.x, goal_.v.y, goal_.v.z));
      Eigen::Vector3d a = enuNedSwap(Eigen::Vector3d(goal_.a.x, goal_.a.y, goal_.a.z));
      Eigen::Vector3d j = enuNedSwap(Eigen::Vector3d(goal_.j.x, goal_.j.y, goal_.j.z));
      msg.velocity = {static_cast<float>(v.x()), static_cast<float>(v.y()),
                      static_cast<float>(v.z())};
      msg.acceleration = {static_cast<float>(a.x()), static_cast<float>(a.y()),
                          static_cast<float>(a.z())};
      msg.jerk = {static_cast<float>(j.x()), static_cast<float>(j.y()), static_cast<float>(j.z())};
      msg.yawspeed = static_cast<float>(-goal_.dyaw);
    }
  } else {
    msg.position = {static_cast<float>(hold_pos_ned_.x()), static_cast<float>(hold_pos_ned_.y()),
                    static_cast<float>(hold_pos_ned_.z())};
    msg.yaw = static_cast<float>(hold_yaw_ned_);
  }

  msg.timestamp = nowMicros();
  pub_trajectory_setpoint_->publish(msg);
}

void Px4OffboardBridge::publishVehicleCommand(uint16_t command, float param1, float param2) {
  VehicleCommand msg{};
  msg.param1 = param1;
  msg.param2 = param2;
  msg.command = command;
  msg.target_system = target_system_;
  msg.target_component = 1;
  msg.source_system = 1;
  msg.source_component = 1;
  msg.from_external = true;
  msg.timestamp = nowMicros();
  pub_vehicle_command_->publish(msg);
}

uint64_t Px4OffboardBridge::nowMicros() { return this->get_clock()->now().nanoseconds() / 1000; }

Eigen::Vector3d Px4OffboardBridge::enuNedSwap(const Eigen::Vector3d& v) {
  return Eigen::Vector3d(v.y(), v.x(), -v.z());
}

double Px4OffboardBridge::yawEnuNedSwap(double yaw) { return wrapPi(M_PI / 2.0 - yaw); }

double Px4OffboardBridge::wrapPi(double x) { return std::atan2(std::sin(x), std::cos(x)); }

Eigen::Vector3d Px4OffboardBridge::globalToLocalNed(double lat, double lon, double alt) const {
  // Azimuthal equidistant projection about the EKF origin (matches PX4's MapProjection::project)
  const double lat_rad = lat * M_PI / 180.0;
  const double lon_rad = lon * M_PI / 180.0;
  const double ref_lat_rad = ref_lat_ * M_PI / 180.0;
  const double ref_lon_rad = ref_lon_ * M_PI / 180.0;

  const double cos_d_lon = std::cos(lon_rad - ref_lon_rad);
  const double arg = std::clamp(std::sin(ref_lat_rad) * std::sin(lat_rad) +
                                    std::cos(ref_lat_rad) * std::cos(lat_rad) * cos_d_lon,
                                -1.0, 1.0);
  const double c = std::acos(arg);
  const double k = (std::fabs(c) > 0.0) ? (c / std::sin(c)) : 1.0;

  const double x = k *
                   (std::cos(ref_lat_rad) * std::sin(lat_rad) -
                    std::sin(ref_lat_rad) * std::cos(lat_rad) * cos_d_lon) *
                   kEarthRadius;
  const double y = k * std::cos(lat_rad) * std::sin(lon_rad - ref_lon_rad) * kEarthRadius;
  const double z = -(alt - ref_alt_);
  return Eigen::Vector3d(x, y, z);
}

int main(int argc, char* argv[]) {
  rclcpp::init(argc, argv);
  auto node = std::make_shared<Px4OffboardBridge>();
  rclcpp::spin(node);
  rclcpp::shutdown();
  return 0;
}
