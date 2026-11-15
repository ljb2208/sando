/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <sando/convert_odom_to_state.hpp>

OdometryToStateNode::OdometryToStateNode() : Node("odometry_to_state_node") {
  // Subscriber for nav_msgs/Odometry
  odometry_sub_ = this->create_subscription<nav_msgs::msg::Odometry>(
      "odom", 10, std::bind(&OdometryToStateNode::callback, this, std::placeholders::_1));

  // Publisher for dynus_interfaces/State
  state_publisher_ = this->create_publisher<dynus_interfaces::msg::State>("state", 10);

  twist_in_body_frame_ = this->declare_parameter("twist_in_body_frame", false);
  publish_tf_ = this->declare_parameter("publish_tf", false);
  restamp_ = this->declare_parameter("restamp", false);
  if (publish_tf_) tf_broadcaster_ = std::make_unique<tf2_ros::TransformBroadcaster>(*this);
}

void OdometryToStateNode::callback(const nav_msgs::msg::Odometry::SharedPtr odom_msg) {
  // Construct the State message
  auto state_msg = dynus_interfaces::msg::State();
  state_msg.header.stamp = restamp_ ? this->now() : rclcpp::Time(odom_msg->header.stamp);
  state_msg.header.frame_id = odom_msg->header.frame_id;

  // Set position from Odometry
  state_msg.pos.x = odom_msg->pose.pose.position.x;
  state_msg.pos.y = odom_msg->pose.pose.position.y;
  state_msg.pos.z = odom_msg->pose.pose.position.z;

  // Set velocity from Odometry (in the pose frame)
  state_msg.vel.x = odom_msg->twist.twist.linear.x;
  state_msg.vel.y = odom_msg->twist.twist.linear.y;
  state_msg.vel.z = odom_msg->twist.twist.linear.z;
  if (twist_in_body_frame_) {
    const auto& q = odom_msg->pose.pose.orientation;
    const Eigen::Quaterniond rot(q.w, q.x, q.y, q.z);
    const Eigen::Vector3d v = rot.normalized() * Eigen::Vector3d(
        state_msg.vel.x, state_msg.vel.y, state_msg.vel.z);
    state_msg.vel.x = v.x();
    state_msg.vel.y = v.y();
    state_msg.vel.z = v.z();
  }

  // Set orientation from Odometry
  state_msg.quat.x = odom_msg->pose.pose.orientation.x;
  state_msg.quat.y = odom_msg->pose.pose.orientation.y;
  state_msg.quat.z = odom_msg->pose.pose.orientation.z;
  state_msg.quat.w = odom_msg->pose.pose.orientation.w;

  // Publish the State message
  state_publisher_->publish(state_msg);

  if (tf_broadcaster_) {
    geometry_msgs::msg::TransformStamped tf;
    tf.header = odom_msg->header;
    tf.header.stamp = state_msg.header.stamp;
    tf.child_frame_id = odom_msg->child_frame_id;
    tf.transform.translation.x = odom_msg->pose.pose.position.x;
    tf.transform.translation.y = odom_msg->pose.pose.position.y;
    tf.transform.translation.z = odom_msg->pose.pose.position.z;
    tf.transform.rotation = odom_msg->pose.pose.orientation;
    tf_broadcaster_->sendTransform(tf);
  }
}

int main(int argc, char* argv[]) {
  rclcpp::init(argc, argv);
  auto node = std::make_shared<OdometryToStateNode>();

  rclcpp::spin(node);
  rclcpp::shutdown();
  return 0;
}
