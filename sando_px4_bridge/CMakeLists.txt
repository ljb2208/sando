cmake_minimum_required(VERSION 3.8)
project(sando_px4_bridge)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(dynus_interfaces REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(px4_msgs REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(include ${EIGEN3_INCLUDE_DIR})

add_executable(px4_offboard_bridge src/px4_offboard_bridge.cpp)
ament_target_dependencies(px4_offboard_bridge rclcpp dynus_interfaces geometry_msgs px4_msgs Eigen3)

# install for executables
install(TARGETS
  px4_offboard_bridge
  DESTINATION lib/${PROJECT_NAME}
)

# install for launch files
install(DIRECTORY
  launch
  DESTINATION share/${PROJECT_NAME}
)

ament_package()
