#!/usr/bin/env python3

# /* ----------------------------------------------------------------------------
#  * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import DeclareLaunchArgument, OpaqueFunction
from launch.substitutions import LaunchConfiguration


def generate_launch_description():
    # Namespace argument
    namespace_arg = DeclareLaunchArgument(
        "namespace", default_value="quadrotor", description="Namespace of the nodes"
    )

    # Parameters for the Px4OffboardBridge node
    param_px4_namespace = DeclareLaunchArgument(
        "px4_namespace",
        default_value="",
        description="PX4 uXRCE-DDS namespace (e.g. 'px4_1' for multi-vehicle), empty for none",
    )
    param_auto_engage = DeclareLaunchArgument(
        "auto_engage",
        default_value="false",
        description="Switch PX4 to offboard and arm automatically once SANDO publishes a goal",
    )
    param_publish_state = DeclareLaunchArgument(
        "publish_state",
        default_value="true",
        description="Publish SANDO state from PX4 vehicle_odometry",
    )

    def launch_setup(context, *args, **kwargs):
        namespace = LaunchConfiguration("namespace").perform(context)
        px4_namespace = LaunchConfiguration("px4_namespace").perform(context)
        px4_prefix = f"/{px4_namespace}/fmu" if px4_namespace else "/fmu"

        # Px4OffboardBridge Node
        node = Node(
            package="sando_px4_bridge",
            executable="px4_offboard_bridge",
            name="px4_offboard_bridge",
            namespace=namespace,
            output="screen",
            parameters=[
                {"auto_engage": LaunchConfiguration("auto_engage")},
                {"publish_state": LaunchConfiguration("publish_state")},
                {"offboard_control_mode_topic": f"{px4_prefix}/in/offboard_control_mode"},
                {"trajectory_setpoint_topic": f"{px4_prefix}/in/trajectory_setpoint"},
                {"vehicle_command_topic": f"{px4_prefix}/in/vehicle_command"},
                {"vehicle_odometry_topic": f"{px4_prefix}/out/vehicle_odometry"},
                {"vehicle_local_position_topic": f"{px4_prefix}/out/vehicle_local_position"},
                {"vehicle_status_topic": f"{px4_prefix}/out/vehicle_status"},
                {"position_setpoint_triplet_topic": f"{px4_prefix}/out/position_setpoint_triplet"},
            ],
            remappings=[
                ("goal", "goal"),  # Incoming SANDO goal (setpoint) topic
                ("state", "state"),  # Outgoing SANDO state topic
                ("term_goal", "term_goal"),  # Outgoing SANDO terminal goal topic
            ],
        )

        return [node]

    return LaunchDescription(
        [
            namespace_arg,
            param_px4_namespace,
            param_auto_engage,
            param_publish_state,
            OpaqueFunction(function=launch_setup),
        ]
    )
