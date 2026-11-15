/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

// Trajectory tracking for an Ackermann-steered rover: follows SANDO's setpoints ("goal": position,
// velocity, acceleration) and publishes the forward speed and the yaw rate the steering produces
// on "cmd_vel" (geometry_msgs/Twist, as Gazebo's AckermannSteering system and most rover bases
// take it), plus the front-wheel angle on "steering_angle" [rad].
//
// Control law (Kanayama-style, in the rover frame):
//   v     = v_ref * cos(e_yaw) + k_x * e_x
//   omega = omega_ref + v_ref * k_y * e_y + |v_ref| * k_yaw * sin(e_yaw)
//   delta = atan(wheelbase * omega / v), clamped to +-max_steering and rate limited
// with e_x, e_y the reference position in the rover frame, e_yaw = yaw_ref - yaw, yaw_ref the
// planned vehicle heading (SANDO's goal yaw), v_ref the reference speed, negative when the
// reference velocity points against the heading (reversing), and omega_ref the goal's dyaw.

#include <algorithm>
#include <cmath>
#include <mutex>

#include <dynus_interfaces/msg/goal.hpp>
#include <dynus_interfaces/msg/state.hpp>
#include <geometry_msgs/msg/twist.hpp>
#include <rclcpp/rclcpp.hpp>
#include <std_msgs/msg/float64.hpp>

using namespace std::chrono_literals;

namespace {
double wrap(double a) { return std::atan2(std::sin(a), std::cos(a)); }

double yawFromQuat(const geometry_msgs::msg::Quaternion& q) {
  return std::atan2(2.0 * (q.w * q.z + q.x * q.y), 1.0 - 2.0 * (q.y * q.y + q.z * q.z));
}
}  // namespace

class AckermannController : public rclcpp::Node {
 public:
  AckermannController() : Node("ackermann_controller") {
    wheelbase_ = declare_parameter("ackermann_wheelbase", 0.5);
    max_steer_ = declare_parameter("ackermann_max_steering_deg", 30.0) * M_PI / 180.0;
    max_steer_rate_ = declare_parameter("max_steering_rate_deg", 120.0) * M_PI / 180.0;
    v_max_ = declare_parameter("v_max", 2.0);
    kx_ = declare_parameter("k_x", 1.0);
    ky_ = declare_parameter("k_y", 1.5);
    kyaw_ = declare_parameter("k_yaw", 2.0);
    rate_hz_ = declare_parameter("rate_hz", 50.0);
    goal_timeout_ = declare_parameter("goal_timeout", 0.5);

    pub_cmd_ = create_publisher<geometry_msgs::msg::Twist>("cmd_vel", 10);
    pub_steer_ = create_publisher<std_msgs::msg::Float64>("steering_angle", 10);
    sub_goal_ = create_subscription<dynus_interfaces::msg::Goal>(
        "goal", 10, [this](const dynus_interfaces::msg::Goal& m) {
          std::lock_guard<std::mutex> lk(mtx_);
          goal_ = m;
          goal_t_ = now();
          have_goal_ = true;
        });
    sub_state_ = create_subscription<dynus_interfaces::msg::State>(
        "state", 10, [this](const dynus_interfaces::msg::State& m) {
          std::lock_guard<std::mutex> lk(mtx_);
          state_ = m;
          have_state_ = true;
        });
    timer_ = create_wall_timer(
        std::chrono::duration<double>(1.0 / rate_hz_), std::bind(&AckermannController::step, this));
    RCLCPP_INFO(
        get_logger(), "Ackermann controller: wheelbase %.2f m, max steering %.1f deg", wheelbase_,
        max_steer_ * 180.0 / M_PI);
  }

 private:
  void step() {
    dynus_interfaces::msg::Goal g;
    dynus_interfaces::msg::State s;
    rclcpp::Time goal_t;
    {
      std::lock_guard<std::mutex> lk(mtx_);
      if (!have_goal_ || !have_state_) return;
      g = goal_;
      s = state_;
      goal_t = goal_t_;
    }
    const double dt = 1.0 / rate_hz_;

    double v_cmd = 0.0, steer = steer_;
    if ((now() - goal_t).seconds() <= goal_timeout_) {  // stop if SANDO stops sending
      const double yaw = yawFromQuat(s.quat);
      const double vx = g.v.x, vy = g.v.y;
      const double speed_ref = std::hypot(vx, vy);
      // SANDO's goal yaw is the planned vehicle heading; the reference velocity points along it
      // when driving forward and against it when reversing
      const double yaw_ref = g.yaw;
      const double dir =
          (speed_ref > 0.05 && vx * std::cos(yaw_ref) + vy * std::sin(yaw_ref) < 0.0) ? -1.0 : 1.0;
      const double v_ref = dir * speed_ref;  // signed
      const double omega_ref = speed_ref > 0.05 ? g.dyaw : 0.0;

      // Reference position in the rover frame
      const double dx = g.p.x - s.pos.x, dy = g.p.y - s.pos.y;
      const double ex = std::cos(yaw) * dx + std::sin(yaw) * dy;
      const double ey = -std::sin(yaw) * dx + std::cos(yaw) * dy;
      const double eyaw = wrap(yaw_ref - yaw);

      // Kanayama-style law, valid in both directions (|v_ref| on the heading term keeps it
      // stabilizing when reversing)
      v_cmd = std::clamp(v_ref * std::cos(eyaw) + kx_ * ex, -v_max_, v_max_);
      if (dir > 0.0) v_cmd = std::max(v_cmd, 0.0);  // no reversing on a forward reference
      if (dir < 0.0) v_cmd = std::min(v_cmd, 0.0);
      const double omega =
          omega_ref + v_ref * ky_ * ey + std::fabs(v_ref) * kyaw_ * std::sin(eyaw);
      if (std::fabs(v_cmd) > 0.05)
        steer = std::atan(wheelbase_ * omega / v_cmd);
      else if (speed_ref > 0.05)
        steer = std::atan(wheelbase_ * omega_ref / v_ref);  // pre-steer while (nearly) stopped
    }

    // Steering limits: angle and rate
    steer = std::clamp(steer, -max_steer_, max_steer_);
    const double max_delta = max_steer_rate_ * dt;
    steer_ += std::clamp(steer - steer_, -max_delta, max_delta);

    geometry_msgs::msg::Twist cmd;
    cmd.linear.x = v_cmd;
    cmd.angular.z = v_cmd * std::tan(steer_) / wheelbase_;
    pub_cmd_->publish(cmd);
    std_msgs::msg::Float64 st;
    st.data = steer_;
    pub_steer_->publish(st);
  }

  double wheelbase_, max_steer_, max_steer_rate_, v_max_, kx_, ky_, kyaw_, rate_hz_, goal_timeout_;
  double steer_ = 0.0;
  std::mutex mtx_;
  dynus_interfaces::msg::Goal goal_;
  dynus_interfaces::msg::State state_;
  rclcpp::Time goal_t_;
  bool have_goal_ = false, have_state_ = false;
  rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr pub_cmd_;
  rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr pub_steer_;
  rclcpp::Subscription<dynus_interfaces::msg::Goal>::SharedPtr sub_goal_;
  rclcpp::Subscription<dynus_interfaces::msg::State>::SharedPtr sub_state_;
  rclcpp::TimerBase::SharedPtr timer_;
};

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<AckermannController>());
  rclcpp::shutdown();
  return 0;
}
