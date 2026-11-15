"""Onboard SANDO launch: planner node, fake sim, and optional obstacle tracker for hardware or simulation."""

import os
from launch import LaunchDescription
from launch_ros.actions import Node
from ament_index_python.packages import get_package_share_directory, PackageNotFoundError
from launch.actions import DeclareLaunchArgument, OpaqueFunction
from launch.substitutions import Command, PathJoinSubstitution, LaunchConfiguration
from launch_ros.substitutions import FindPackageShare
from launch_ros.parameter_descriptions import ParameterValue
import yaml
from math import radians


def gazebo_backend():
    """Gazebo Classic ("classic") when gazebo_ros is installed (Humble), else Gazebo Sim ("gz")."""
    try:
        get_package_share_directory("gazebo_ros")
        return "classic"
    except PackageNotFoundError:
        return "gz"


def convert_str_to_bool(str):
    return True if (str == "true" or str == "True" or str == 1 or str == "1") else False


def generate_launch_description():

    # Declare launch arguments

    # initial position and yaw of the quadrotor
    x_arg = DeclareLaunchArgument(
        "x", default_value="20.0", description="Initial x position of the quadrotor"
    )
    y_arg = DeclareLaunchArgument(
        "y", default_value="9.0", description="Initial y position of the quadrotor"
    )
    z_arg = DeclareLaunchArgument(
        "z", default_value="2.0", description="Initial z position of the quadrotor"
    )
    yaw_arg = DeclareLaunchArgument(
        "yaw", default_value="180", description="Initial yaw angle of the quadrotor"
    )
    namespace_arg = DeclareLaunchArgument(
        "namespace", default_value="NX01", description="Namespace of the nodes"
    )  # namespace
    use_obstacle_tracker_arg = DeclareLaunchArgument(
        "use_obstacle_tracker",
        default_value="false",
        description="Flag to indicate whether to start the obstacle tracker node",
    )  # flag to indicate whether to start the obstacle tracker node
    data_file_arg = DeclareLaunchArgument(
        "data_file",
        default_value="",
        description="File name to store data",
    )  # file name to store data
    global_planner_arg = DeclareLaunchArgument(
        "global_planner", default_value="sjps", description="Global planner to use"
    )  # global planner
    solver_arg = DeclareLaunchArgument(
        "solver",
        default_value="",
        description="Trajectory optimization backend: highs_shot, gurobi or auto "
        "(empty: the 'solver' parameter of the config file)",
    )  # trajectory solver
    vehicle_arg = DeclareLaunchArgument(
        "vehicle",
        default_value="uav",
        description="uav, or ackermann (rover: applies config/sando_ackermann.yaml, starts "
        "ackermann_controller and simulates a kinematic bicycle in fake_sim)",
    )
    web_arg = DeclareLaunchArgument(
        "web",
        default_value="false",
        description="Start the SANDO web interface (status, performance, parameters, 3D view)",
    )
    web_port_arg = DeclareLaunchArgument(
        "web_port", default_value="8080", description="Port of the web interface"
    )
    use_benchmark_arg = DeclareLaunchArgument(
        "use_benchmark",
        default_value="false",
        description="Flag to indicate whether to use the global planner benchmark",
    )  # global planner benchmark
    use_hardware_arg = DeclareLaunchArgument(
        "use_hardware",
        default_value="false",
        description="Flag to indicate whether to use hardware or simulation",
    )  # flag to indicte if this is hardware or simulation
    use_onboard_localization_arg = DeclareLaunchArgument(
        "use_onboard_localization",
        default_value="false",
        description="Use onboard localization (DLIO odom) instead of Vicon",
    )
    sim_env_arg = DeclareLaunchArgument(
        "sim_env",
        default_value="",
        description="Simulation environment (gazebo, fake_sim). Empty string uses value from config file.",
    )  # override sim_env from config
    environment_assumption_arg = DeclareLaunchArgument(
        "environment_assumption",
        default_value="",
        description="Override environment_assumption from config (static, dynamic, dynamic_worst_case). Empty string uses value from config file.",
    )
    publish_odom_arg = DeclareLaunchArgument("publish_odom", default_value="true")
    odom_topic_arg = DeclareLaunchArgument(
        "odom_topic", default_value="visual_slam/odom"
    )
    odom_frame_id_arg = DeclareLaunchArgument("odom_frame_id", default_value="map")
    v_max_arg = DeclareLaunchArgument(
        "v_max",
        default_value="",
        description="Override v_max from config. Empty string uses value from config file.",
    )
    hover_avoidance_enabled_arg = DeclareLaunchArgument(
        "hover_avoidance_enabled",
        default_value="",
        description="Override hover_avoidance_enabled from config (true/false). Empty string uses value from config file.",
    )
    ignore_other_trajs_arg = DeclareLaunchArgument(
        "ignore_other_trajs",
        default_value="",
        description="Override ignore_other_trajs from config (true/false). Empty string uses value from config file.",
    )
    skip_initial_yawing_arg = DeclareLaunchArgument(
        "skip_initial_yawing",
        default_value="",
        description="Skip YAWING state and go directly to TRAVELING (true/false). Empty string uses value from config file.",
    )

    # Need to be the same as simulartor.launch.py
    map_size_x_arg = DeclareLaunchArgument("map_size_x", default_value="20.0")
    map_size_y_arg = DeclareLaunchArgument("map_size_y", default_value="20.0")
    map_size_z_arg = DeclareLaunchArgument("map_size_z", default_value="6.0")
    odometry_topic_arg = DeclareLaunchArgument(
        "odometry_topic", default_value="visual_slam/odom"
    )

    # Opaque function to launch nodes
    def launch_setup(context, *args, **kwargs):
        backend = gazebo_backend()

        x = LaunchConfiguration("x").perform(context)
        y = LaunchConfiguration("y").perform(context)
        z = LaunchConfiguration("z").perform(context)
        yaw = LaunchConfiguration("yaw").perform(context)
        namespace = LaunchConfiguration("namespace").perform(context)
        use_obstacle_tracker = convert_str_to_bool(
            LaunchConfiguration("use_obstacle_tracker").perform(context)
        )
        data_file = LaunchConfiguration("data_file").perform(context)
        global_planner = LaunchConfiguration("global_planner").perform(context)
        solver = LaunchConfiguration("solver").perform(context)
        web = convert_str_to_bool(LaunchConfiguration("web").perform(context))
        vehicle = LaunchConfiguration("vehicle").perform(context)
        if vehicle not in ("uav", "ackermann"):
            raise RuntimeError(f"vehicle:={vehicle}: use uav or ackermann")
        web_port = int(LaunchConfiguration("web_port").perform(context))
        use_benchmark = convert_str_to_bool(
            LaunchConfiguration("use_benchmark").perform(context)
        )
        use_hardware = convert_str_to_bool(
            LaunchConfiguration("use_hardware").perform(context)
        )
        use_onboard_localization = convert_str_to_bool(
            LaunchConfiguration("use_onboard_localization").perform(context)
        )
        sim_env_override = LaunchConfiguration("sim_env").perform(context)
        environment_assumption_override = LaunchConfiguration(
            "environment_assumption"
        ).perform(context)
        publish_odom = convert_str_to_bool(
            LaunchConfiguration("publish_odom").perform(context)
        )
        odom_topic = LaunchConfiguration("odom_topic").perform(context)
        odom_frame_id = LaunchConfiguration("odom_frame_id").perform(context)
        v_max_override = LaunchConfiguration("v_max").perform(context)
        hover_avoidance_enabled_override = LaunchConfiguration(
            "hover_avoidance_enabled"
        ).perform(context)
        ignore_other_trajs_override = LaunchConfiguration("ignore_other_trajs").perform(
            context
        )
        skip_initial_yawing_override = LaunchConfiguration(
            "skip_initial_yawing"
        ).perform(context)
        base_frame_id = namespace + "/base_link"
        map_size_x = float(LaunchConfiguration("map_size_x").perform(context))
        map_size_y = float(LaunchConfiguration("map_size_y").perform(context))
        map_size_z = float(LaunchConfiguration("map_size_z").perform(context))
        odometry_topic = LaunchConfiguration("odometry_topic").perform(context)

        # The path to the urdf file
        urdf_path = PathJoinSubstitution(
            [FindPackageShare("sando"), "urdf", "quadrotor.urdf.xacro"]
        )
        if use_hardware:
            parameters_path = os.path.join(
                get_package_share_directory("sando"),
                "config",
                "sando_hw_quadrotor.yaml",
            )
        else:
            parameters_path = os.path.join(
                get_package_share_directory("sando"), "config", "sando.yaml"
            )

        # Get the dict of parameters from the yaml file
        with open(parameters_path, "r") as file:
            parameters = yaml.safe_load(file)

        # Extract specific node parameters
        parameters = parameters["sando_node"]["ros__parameters"]

        # Rover: overlay the Ackermann parameters and start on the planning plane
        ackermann = vehicle == "ackermann"
        if ackermann:
            with open(
                os.path.join(
                    get_package_share_directory("sando"), "config", "sando_ackermann.yaml"
                ),
                "r",
            ) as file:
                parameters.update(yaml.safe_load(file)["sando_node"]["ros__parameters"])
            urdf_path = PathJoinSubstitution(
                [FindPackageShare("sando"), "urdf", "rover.urdf.xacro"]
            )

        # Override sim_env if provided (hardware configs may not have sim_env)
        if sim_env_override:
            parameters["sim_env"] = sim_env_override
        elif "sim_env" not in parameters:
            parameters["sim_env"] = ""

        # Override environment_assumption if provided
        if environment_assumption_override:
            parameters["environment_assumption"] = environment_assumption_override

        # Override v_max if provided
        if v_max_override:
            parameters["v_max"] = float(v_max_override)

        # Override hover_avoidance_enabled if provided
        if hover_avoidance_enabled_override:
            parameters["hover_avoidance_enabled"] = convert_str_to_bool(
                hover_avoidance_enabled_override
            )

        # Override ignore_other_trajs if provided
        if ignore_other_trajs_override:
            parameters["ignore_other_trajs"] = convert_str_to_bool(
                ignore_other_trajs_override
            )

        # Override skip_initial_yawing if provided
        if skip_initial_yawing_override:
            parameters["skip_initial_yawing"] = convert_str_to_bool(
                skip_initial_yawing_override
            )

        if solver:
            parameters["solver"] = solver

        # (after the sim_env override) In Gazebo the rover drives on the ground with real physics;
        # in fake_sim it starts on the planning plane
        spawn_z = z
        if ackermann:
            z = str(parameters["ground_plane_z"])
            spawn_z = "0.25" if parameters["sim_env"] == "gazebo" else z
        gazebo_rover = ackermann and parameters["sim_env"] == "gazebo"

        # Update parameters for benchmarking
        parameters["file_path"] = data_file
        parameters["use_benchmark"] = bool(use_benchmark)
        if use_benchmark:
            parameters["global_planner"] = global_planner

        # Create a Sando node
        sando_node = Node(
            package="sando",
            executable="sando",
            name="sando_node",
            namespace=namespace,
            output="screen",
            emulate_tty=True,
            parameters=[parameters],
            # prefix='xterm -e gdb -q -ex run --args', # gdb debugging
            arguments=["--ros-args", "--log-level", "error"],
        )

        # Robot state publisher node
        robot_state_publisher_node = Node(
            package="robot_state_publisher",
            executable="robot_state_publisher",
            name="robot_state_publisher",
            output="screen",
            namespace=namespace,
            parameters=[
                {
                    "robot_description": ParameterValue(
                        Command(
                            [
                                "xacro ",
                                urdf_path,
                                " namespace:=",
                                namespace,
                                " d435_range_max_depth:=",
                                str(parameters["depth_camera_depth_max"]),
                                " gazebo:=",
                                backend,
                            ]
                        ),
                        value_type=str,
                    ),
                    "use_sim_time": False,
                    "frame_prefix": namespace + "/",
                }
            ],
            arguments=["--ros-args", "--log-level", "error"],
        )

        # Spawn entity node for Gazebo
        # Get the start position and yaw from the parameters
        yaw = str(radians(float(yaw)))
        if backend == "gz":
            spawn_entity_node = Node(
                package="ros_gz_sim",
                executable="create",
                name="spawn_entity",
                namespace=namespace,
                arguments=[
                    "-world",
                    "default",
                    "-topic",
                    "robot_description",
                    "-name",
                    namespace,
                    "-x",
                    x,
                    "-y",
                    y,
                    "-z",
                    spawn_z,
                    "-Y",
                    yaw,
                ],
            )

            # Bridge the D435 and Mid-360 topics from Gazebo Sim to the names the Classic plugin used
            d435 = f"/{namespace}/d435"
            gz_bridge_node = Node(
                package="ros_gz_bridge",
                executable="parameter_bridge",
                name="gz_sensor_bridge",
                namespace=namespace,
                output="log",
                arguments=[
                    f"{d435}/depth/color/points@sensor_msgs/msg/PointCloud2[gz.msgs.PointCloudPacked",
                    f"{d435}/depth/color@sensor_msgs/msg/Image[gz.msgs.Image",
                    f"{d435}/depth/camera_info@sensor_msgs/msg/CameraInfo[gz.msgs.CameraInfo",
                    f"{d435}/color/image_raw@sensor_msgs/msg/Image[gz.msgs.Image",
                    f"{d435}/color/camera_info@sensor_msgs/msg/CameraInfo[gz.msgs.CameraInfo",
                    f"/{namespace}/mid360/points@sensor_msgs/msg/PointCloud2[gz.msgs.PointCloudPacked",
                ]
                + (
                    [
                        # Rover: commands in, ground-truth odometry and wheel joints out
                        f"/{namespace}/cmd_vel@geometry_msgs/msg/Twist]gz.msgs.Twist",
                        f"/{namespace}/odometry_gt@nav_msgs/msg/Odometry[gz.msgs.Odometry",
                        f"/{namespace}/joint_states@sensor_msgs/msg/JointState[gz.msgs.Model",
                    ]
                    if gazebo_rover
                    else []
                ),
                remappings=[
                    (f"{d435}/depth/color", f"{d435}/depth/image_raw"),
                    (f"/{namespace}/mid360/points", f"/{namespace}/mid360_PointCloud2"),
                ],
            )
        else:
            spawn_entity_node = Node(
                package="gazebo_ros",
                executable="spawn_entity.py",
                name="spawn_entity",
                namespace=namespace,
                parameters=[
                    {
                        "use_sim_time": False,
                    }
                ],
                arguments=[
                    "-topic",
                    "robot_description",
                    "-entity",
                    namespace,
                    "-x",
                    x,
                    "-y",
                    y,
                    "-z",
                    z,
                    "-Y",
                    yaw,
                ],
            )

        # Create an obstacle tracker node
        obstacle_tracker_node = Node(
            package="sando",
            executable="obstacle_tracker_node",
            namespace=namespace,
            name="obstacle_tracker_node",
            emulate_tty=True,
            parameters=[parameters],
            # prefix='xterm -e gdb -ex run --args', # gdb debugging
            output="screen",
            remappings=[("point_cloud", "d435/depth/color/points")],
        )

        # Convert pose and twist (from Vicon) to state
        pose_twist_to_state_node = Node(
            package="sando",
            executable="convert_vicon_to_state",
            name="convert_vicon_to_state",
            namespace=namespace,
            remappings=[
                ("world", "world"),  # Remap incoming PoseStamped topic
                ("twist", "twist"),  # Remap incoming TwistStamped topic
                ("state", "state"),  # Remap outgoing State topic
            ],
            emulate_tty=True,
            output="screen",
            # prefix='xterm -e gdb -ex run --args', # gdb debugging
            # arguments=['--ros-args', '--log-level', 'error']
        )

        # Convert odom (from DLIO) to global-frame state and pose via TF2
        odom_to_global_state_node = Node(
            package="sando",
            executable="odom_to_global_state",
            name="odom_to_global_state",
            namespace=namespace,
            remappings=[
                ("odom", "dlio/odom_node/odom"),
                ("state", "state"),
                ("global_pose", "global_pose"),
            ],
            emulate_tty=True,
            output="screen",
        )

        # When using ground robot, we don't need to send the exact state to gazebo - the state will be taken care of by wheel controllers
        # send_state_to_gazebo = False if use_ground_robot else True
        # Create a fake sim node
        fake_sim_node = Node(
            package="sando",
            executable="fake_sim",
            name="fake_sim",
            namespace=namespace,
            emulate_tty=True,
            parameters=[
                {
                    "start_pos": [float(x), float(y), float(z)],
                    "start_yaw": float(yaw),
                    "send_state_to_gazebo": parameters["sim_env"] == "gazebo",
                    "gazebo_backend": backend,
                    "visual_level": parameters["visual_level"],
                    "publish_odom": publish_odom,
                    "odom_topic": odom_topic,
                    "odom_frame_id": odom_frame_id,
                    "base_frame_id": base_frame_id,
                    "vehicle_model": "ackermann" if ackermann else "ideal",
                    "ackermann_wheelbase": float(parameters.get("ackermann_wheelbase", 0.5)),
                    "ackermann_max_steering_deg": float(
                        parameters.get("ackermann_max_steering_deg", 30.0)
                    ),
                    "ackermann_max_accel": float(parameters.get("a_max", 2.0)),
                }
            ],
            # prefix='xterm -e gdb -q -ex run --args', # gdb debugging
            output="screen",
        )

        camera_file = os.path.join(
            get_package_share_directory("local_sensing"), "config", "camera.yaml"
        )

        pcl_render_node = Node(
            package="local_sensing",
            executable="pcl_render_node",
            namespace=namespace,
            name="pcl_render_node",
            output="screen",
            parameters=[
                {"sensing_horizon": 5.0},
                {"sensing_rate": 30.0},
                {"estimation_rate": 30.0},
                {"map/x_size": map_size_x},
                {"map/y_size": map_size_y},
                {"map/z_size": map_size_z},
                camera_file,
            ],
            remappings=[
                ("global_map", "/map_generator/global_cloud"),
                ("odometry", odometry_topic),
                ("depth", "pcl_render_node/depth"),
            ],
            # prefix='xterm -e gdb -q -ex run --args', # gdb debugging
        )

        # Return launch description
        nodes_to_start = [sando_node]
        if use_hardware and use_onboard_localization:
            nodes_to_start.append(odom_to_global_state_node)
        elif use_hardware:
            nodes_to_start.append(pose_twist_to_state_node)
        elif gazebo_rover:
            # Rover in Gazebo: physics drives it; state and the map -> base_link TF come from the
            # ground-truth odometry (body-frame twist, sim-time stamps)
            nodes_to_start += [
                robot_state_publisher_node,
                spawn_entity_node,
                gz_bridge_node,
                Node(
                    package="sando",
                    executable="convert_odom_to_state",
                    name="odom_to_state",
                    namespace=namespace,
                    output="screen",
                    parameters=[
                        {"twist_in_body_frame": True, "publish_tf": True, "restamp": True}
                    ],
                    remappings=[("odom", "odometry_gt")],
                ),
            ]
        else:
            nodes_to_start.append(fake_sim_node)
            if parameters["sim_env"] == "gazebo":
                nodes_to_start.append(robot_state_publisher_node)
                nodes_to_start.append(spawn_entity_node)
                if backend == "gz":
                    nodes_to_start.append(gz_bridge_node)
            elif parameters["sim_env"] == "fake_sim":
                nodes_to_start.append(pcl_render_node)
        nodes_to_start.append(obstacle_tracker_node) if use_obstacle_tracker else None
        if ackermann:
            # Tracks SANDO's setpoints: goal + state -> cmd_vel (speed, yaw rate)
            nodes_to_start.append(
                Node(
                    package="sando",
                    executable="ackermann_controller",
                    name="ackermann_controller",
                    namespace=namespace,
                    output="screen",
                    parameters=[
                        {
                            "ackermann_wheelbase": float(parameters["ackermann_wheelbase"]),
                            "ackermann_max_steering_deg": float(
                                parameters["ackermann_max_steering_deg"]
                            ),
                            "v_max": float(parameters["v_max"]),
                        }
                    ],
                )
            )
        if web:
            nodes_to_start.append(
                Node(
                    package="sando",
                    executable="sando_web.py",
                    name="sando_web",
                    namespace=namespace,
                    output="screen",
                    parameters=[
                        {
                            "namespace": namespace,
                            "port": web_port,
                            "default_goal_z": float(parameters.get("default_goal_z", 3.0)),
                        }
                    ],
                )
            )

        return nodes_to_start

    # Create launch description
    return LaunchDescription(
        [
            x_arg,
            y_arg,
            z_arg,
            yaw_arg,
            namespace_arg,
            use_obstacle_tracker_arg,
            data_file_arg,
            global_planner_arg,
            solver_arg,
            vehicle_arg,
            web_arg,
            web_port_arg,
            use_benchmark_arg,
            use_hardware_arg,
            use_onboard_localization_arg,
            sim_env_arg,
            environment_assumption_arg,
            publish_odom_arg,
            odom_topic_arg,
            odom_frame_id_arg,
            map_size_x_arg,
            map_size_y_arg,
            map_size_z_arg,
            odometry_topic_arg,
            v_max_arg,
            hover_avoidance_enabled_arg,
            ignore_other_trajs_arg,
            skip_initial_yawing_arg,
            OpaqueFunction(function=launch_setup),
        ]
    )
