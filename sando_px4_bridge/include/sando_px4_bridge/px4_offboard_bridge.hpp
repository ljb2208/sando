/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <optional>
#include "rclcpp/rclcpp.hpp"
#include "dynus_interfaces/msg/goal.hpp"
#include "dynus_interfaces/msg/state.hpp"
#include "geometry_msgs/msg/pose_stamped.hpp"
#include "px4_msgs/msg/offboard_control_mode.hpp"
#include "px4_msgs/msg/position_setpoint_triplet.hpp"
#include "px4_msgs/msg/trajectory_setpoint.hpp"
#include "px4_msgs/msg/vehicle_command.hpp"
#include "px4_msgs/msg/vehicle_local_position.hpp"
#include "px4_msgs/msg/vehicle_odometry.hpp"
#include "px4_msgs/msg/vehicle_status.hpp"

/** @brief ROS 2 node bridging SANDO and PX4 over the uXRCE-DDS interface.
 *
 * - Reads the navigator's current position setpoint (mission waypoint / QGC "Go to") from PX4
 *   and forwards it to SANDO as the terminal goal.
 * - Converts PX4 odometry (NED/FRD) into SANDO's State (ENU/FLU).
 * - Streams SANDO's Goal back to PX4 as OffboardControlMode + TrajectorySetpoint (NED).
 */
class Px4OffboardBridge : public rclcpp::Node {
 public:
  /** @brief Construct the node, load parameters, and set up ROS interfaces. */
  Px4OffboardBridge();

 private:
  void odometryCallback(const px4_msgs::msg::VehicleOdometry::SharedPtr msg);
  void localPositionCallback(const px4_msgs::msg::VehicleLocalPosition::SharedPtr msg);
  void vehicleStatusCallback(const px4_msgs::msg::VehicleStatus::SharedPtr msg);
  void setpointTripletCallback(const px4_msgs::msg::PositionSetpointTriplet::SharedPtr msg);
  void goalCallback(const dynus_interfaces::msg::Goal::SharedPtr msg);
  void offboardTimerCallback();

  void publishOffboardControlMode();
  void publishTrajectorySetpoint();
  void publishVehicleCommand(uint16_t command, float param1 = 0.0, float param2 = 0.0);
  uint64_t nowMicros();

  // Frame conversions (ENU <-> NED is its own inverse)
  static Eigen::Vector3d enuNedSwap(const Eigen::Vector3d& v);
  static double yawEnuNedSwap(double yaw);
  static double wrapPi(double x);
  Eigen::Vector3d globalToLocalNed(double lat, double lon, double alt) const;

  // Publishers and Subscribers (PX4 side)
  rclcpp::Publisher<px4_msgs::msg::OffboardControlMode>::SharedPtr pub_offboard_control_mode_;
  rclcpp::Publisher<px4_msgs::msg::TrajectorySetpoint>::SharedPtr pub_trajectory_setpoint_;
  rclcpp::Publisher<px4_msgs::msg::VehicleCommand>::SharedPtr pub_vehicle_command_;
  rclcpp::Subscription<px4_msgs::msg::VehicleOdometry>::SharedPtr sub_odometry_;
  rclcpp::Subscription<px4_msgs::msg::VehicleLocalPosition>::SharedPtr sub_local_position_;
  rclcpp::Subscription<px4_msgs::msg::VehicleStatus>::SharedPtr sub_vehicle_status_;
  rclcpp::Subscription<px4_msgs::msg::PositionSetpointTriplet>::SharedPtr sub_setpoint_triplet_;

  // Publishers and Subscribers (SANDO side)
  rclcpp::Publisher<dynus_interfaces::msg::State>::SharedPtr pub_state_;
  rclcpp::Publisher<geometry_msgs::msg::PoseStamped>::SharedPtr pub_term_goal_;
  rclcpp::Subscription<dynus_interfaces::msg::Goal>::SharedPtr sub_goal_;

  // Timers
  rclcpp::TimerBase::SharedPtr timer_;

  // Parameters
  std::string frame_id_;
  bool publish_state_;
  bool auto_engage_;
  double setpoint_rate_hz_;
  double goal_timeout_;
  double term_goal_update_threshold_;
  uint8_t target_system_;

  // Vehicle state (NED)
  Eigen::Vector3d pos_ned_;
  double yaw_ned_;
  uint8_t nav_state_;
  uint8_t arming_state_;

  // Local position reference (for lat/lon/alt -> local NED)
  double ref_lat_;
  double ref_lon_;
  double ref_alt_;

  // Latest SANDO goal and the hover setpoint used when no fresh goal is available (NED)
  dynus_interfaces::msg::Goal goal_;
  rclcpp::Time goal_stamp_;
  Eigen::Vector3d hold_pos_ned_;
  double hold_yaw_ned_;
  std::optional<Eigen::Vector3d> last_term_goal_ned_;

  // Counters
  int setpoint_counter_;

  // Flags
  bool odom_initialized_;
  bool ref_initialized_;
  bool goal_initialized_;
  bool engage_sent_;
};
