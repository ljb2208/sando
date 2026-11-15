/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#pragma once

#include <dynus_interfaces/msg/state.hpp>
#include "sando/sando.hpp"
#include "rclcpp/rclcpp.hpp"
#include "geometry_msgs/msg/quaternion.hpp"
#include "geometry_msgs/msg/vector3.hpp"
#include "nav_msgs/msg/odometry.hpp"
#include "std_msgs/msg/header.hpp"
#include "tf2_ros/transform_broadcaster.h"

/** @brief ROS 2 node that converts nav_msgs/Odometry messages into dynus_interfaces/State messages.
 */
class OdometryToStateNode : public rclcpp::Node {
 public:
  /** @brief Construct the node and set up the odometry subscriber and state publisher. */
  OdometryToStateNode();

 private:
  void callback(const nav_msgs::msg::Odometry::SharedPtr odom_msg);

  rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odometry_sub_;
  rclcpp::Publisher<dynus_interfaces::msg::State>::SharedPtr state_publisher_;

  // twist_in_body_frame: the odometry twist is in the child (body) frame, as nav_msgs specifies
  // and Gazebo publishes it; rotate it into the pose frame for the state
  bool twist_in_body_frame_ = false;
  // publish_tf: broadcast pose frame -> child frame (e.g. map -> NX01/base_link for the mapper)
  bool publish_tf_ = false;
  // restamp: stamp with this node's clock (e.g. Gazebo sim-time odometry in a wall-clock system)
  bool restamp_ = false;
  std::unique_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_;
};