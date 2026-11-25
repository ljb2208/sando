#!/bin/bash
# --------------------------------------------------------------------------
# SANDO Build Script
# Builds only the ROS 2 (colcon) components: DecompROS2 in its own workspace, then the SANDO
# workspace (sando, the deps/ packages and sando_px4_bridge). Installs nothing system-wide; run
# setup.sh first on a new machine for ROS 2, Gurobi/HiGHS/SHOT, Livox-SDK2, etc.
#
# Usage:
#   ./build.sh [options]
#
#   -d, --distro NAME     ROS 2 distro in /opt/ros (default: newest installed)
#   -j, --jobs N          Parallel compile jobs (default: 4)
#   -m, --mem-limit SIZE  Memory cap for the build, e.g. 9G; "none" to disable (default: 9G).
#                         Full parallel builds can exhaust RAM and get the desktop session
#                         (e.g. VS Code) killed by systemd-oomd
#   -p, --packages "A B"  Build only these packages of the SANDO workspace (skips DecompROS2)
#       --px4-ws PATH     Workspace providing px4_msgs for sando_px4_bridge
#                         (default: ~/ros2_px4_ws; the bridge is skipped if px4_msgs is missing)
#       --no-gurobi       Do not build the Gurobi solver backend (built when Gurobi is found)
#       --no-highs-shot   Do not build the HiGHS + SHOT solver backend (built when found)
#                         The backend is chosen at launch: solver:=highs_shot|gurobi|auto
#       --clean           Remove build/, install/ and log/ of both workspaces first
#   -h, --help            Show this help
# --------------------------------------------------------------------------

set -e

usage() { sed -n '3,/^# ----*$/p' "$0" | sed '$d; s/^# \{0,1\}//'; }

# Parse arguments
ROS_DISTRO_ARG=""
NUM_JOBS=4
MEM_LIMIT="9G"
PACKAGES=""
PX4_WS="$HOME/ros2_px4_ws"
WITH_GUROBI=ON
WITH_HIGHS_SHOT=ON
CLEAN=false
while [[ $# -gt 0 ]]; do
    case $1 in
        -d|--distro) ROS_DISTRO_ARG="$2"; shift 2 ;;
        -j|--jobs) NUM_JOBS="$2"; shift 2 ;;
        -m|--mem-limit) MEM_LIMIT="$2"; shift 2 ;;
        -p|--packages) PACKAGES="$2"; shift 2 ;;
        --px4-ws) PX4_WS="$2"; shift 2 ;;
        --no-gurobi) WITH_GUROBI=OFF; shift ;;
        --no-highs-shot) WITH_HIGHS_SHOT=OFF; shift ;;
        --clean) CLEAN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SANDO_WS="$(cd "$SCRIPT_DIR/../.." && pwd)"
CODE_DIR="$(cd "$SANDO_WS/.." && pwd)"
DECOMP_WS="$CODE_DIR/decomp_ws"

# Re-run inside a memory-capped systemd scope, so running out of memory only kills the build
if [ "$MEM_LIMIT" != "none" ] && [ -z "$SANDO_BUILD_SCOPED" ] && command -v systemd-run > /dev/null; then
    export SANDO_BUILD_SCOPED="$MEM_LIMIT"
    FLAGS=()
    [ "$WITH_GUROBI" = OFF ] && FLAGS+=(--no-gurobi)
    [ "$WITH_HIGHS_SHOT" = OFF ] && FLAGS+=(--no-highs-shot)
    $CLEAN && FLAGS+=(--clean)
    exec systemd-run --user --scope --quiet -p MemoryMax="$MEM_LIMIT" -p MemorySwapMax=4G \
        "$0" -d "$ROS_DISTRO_ARG" -j "$NUM_JOBS" -m none -p "$PACKAGES" --px4-ws "$PX4_WS" \
        "${FLAGS[@]}"
fi

# Pick the ROS 2 distro
if [ -z "$ROS_DISTRO_ARG" ]; then
    ROS_DISTRO_ARG=$(ls /opt/ros 2>/dev/null | sort | tail -1)
fi
if [ ! -f "/opt/ros/$ROS_DISTRO_ARG/setup.bash" ]; then
    echo "ROS 2 distro '$ROS_DISTRO_ARG' not found in /opt/ros"
    exit 1
fi

# Clean environment: only /opt/ros/<distro> plus the workspaces built here. Conda's python breaks
# ament (no catkin_pkg/em), and overlays sourced from .bashrc (another distro or workspace) get
# mixed into the build
PATH=$(echo "$PATH" | tr ':' '\n' | grep -v -E "conda|/opt/ros/" | paste -sd:)
unset CONDA_PREFIX CONDA_DEFAULT_ENV PYTHONPATH AMENT_PREFIX_PATH CMAKE_PREFIX_PATH \
      COLCON_PREFIX_PATH LD_LIBRARY_PATH ROS_DISTRO ROS_VERSION ROS_PYTHON_VERSION
export PATH
# shellcheck disable=SC1090
source "/opt/ros/$ROS_DISTRO_ARG/setup.bash"

# Each compile job gets its own make -j, so also cap make by load to avoid jobs x jobs processes
export MAKEFLAGS="-j$NUM_JOBS -l$NUM_JOBS"
CMAKE_ARGS=(-DCMAKE_BUILD_TYPE=Release
            -DCMAKE_CXX_COMPILER=/usr/bin/g++
            -DCMAKE_C_COMPILER=/usr/bin/gcc
            -DPython3_EXECUTABLE=/usr/bin/python3)

echo "============================================="
echo "  SANDO Build (ROS 2 components)"
echo "  Distro:    $ROS_DISTRO_ARG"
echo "  Workspace: $SANDO_WS"
echo "  Decomp:    $DECOMP_WS"
echo "  Jobs:      $NUM_JOBS   Memory cap: ${SANDO_BUILD_SCOPED:-none}"
echo "============================================="

if $CLEAN; then
    echo "Cleaning build/, install/ and log/..."
    rm -rf "$DECOMP_WS"/{build,install,log} "$SANDO_WS"/{build,install,log}
fi

# A workspace built for another distro can't be overlaid: refuse instead of mixing them
for ws in "$DECOMP_WS" "$SANDO_WS"; do
    stamp="$ws/install/.ros_distro"
    if [ -f "$stamp" ] && [ "$(cat "$stamp")" != "$ROS_DISTRO_ARG" ]; then
        echo "$ws was built for ROS 2 $(cat "$stamp"), not $ROS_DISTRO_ARG. Re-run with --clean."
        exit 1
    fi
done

# ============================================================
# 1. Workspace symlinks for colcon discovery (same as setup.sh step 7)
# ============================================================
echo ""
echo "STEP 1: Workspace symlinks"
WS_SRC="$SANDO_WS/src"
for d in "$SCRIPT_DIR"/deps/*/; do
    d="${d%/}"
    name=$(basename "$d")
    # Skip deps built in separate workspaces or with plain cmake
    case "$name" in DecompROS2) continue;; esac
    # Gazebo Classic is only available up to Humble: keep just gazebo_msgs and skip the
    # Gazebo Classic simulation packages
    if [ "$ROS_DISTRO_ARG" != "humble" ]; then
        case "$name" in
            gazebo_ros_pkgs)
                [ -e "$WS_SRC/gazebo_msgs" ] || ln -s "$d/gazebo_msgs" "$WS_SRC/gazebo_msgs"
                continue;;
            livox_laser_simulation_ros2|realsense_gazebo_plugin|unitree-go2-ros2) continue;;
        esac
    fi
    if [ ! -e "$WS_SRC/$name" ]; then
        ln -s "$d" "$WS_SRC/$name"
        echo "  Linked: $name"
    fi
done
# sando_px4_bridge lives inside the sando package, which colcon does not search below
[ -e "$WS_SRC/sando_px4_bridge" ] || ln -s "$SCRIPT_DIR/sando_px4_bridge" "$WS_SRC/sando_px4_bridge"

# ============================================================
# 2. DecompROS2 (separate workspace, as in setup.sh)
# ============================================================
if [ -z "$PACKAGES" ]; then
    echo ""
    echo "STEP 2: DecompROS2"
    mkdir -p "$DECOMP_WS/src"
    [ -e "$DECOMP_WS/src/DecompROS2" ] || ln -s "$SCRIPT_DIR/deps/DecompROS2" "$DECOMP_WS/src/DecompROS2"
    cd "$DECOMP_WS"
    # decomp_util first: the other decomp packages use it without declaring the dependency
    colcon build --packages-select decomp_util --parallel-workers "$NUM_JOBS" \
        --cmake-args "${CMAKE_ARGS[@]}"
    source "$DECOMP_WS/install/local_setup.bash"
    colcon build --parallel-workers "$NUM_JOBS" --cmake-args "${CMAKE_ARGS[@]}"
    echo "$ROS_DISTRO_ARG" > "$DECOMP_WS/install/.ros_distro"
else
    echo ""
    echo "STEP 2: DecompROS2 (skipped: --packages given)"
fi
if [ ! -f "$DECOMP_WS/install/local_setup.bash" ]; then
    echo "DecompROS2 is not built yet: run ./build.sh without --packages first"
    exit 1
fi
source "$DECOMP_WS/install/local_setup.bash"

# ============================================================
# 3. SANDO workspace
# ============================================================
echo ""
echo "STEP 3: SANDO workspace"
SKIP_ARGS=()
# px4_msgs for sando_px4_bridge (local_setup: setup.bash of an overlay re-sources its underlay)
if [ -f "$PX4_WS/install/local_setup.bash" ]; then
    source "$PX4_WS/install/local_setup.bash"
fi
PX4_CACHE="$PX4_WS/build/px4_msgs/CMakeCache.txt"
if ! ros2 pkg prefix px4_msgs > /dev/null 2>&1; then
    echo "  px4_msgs not found (looked in $PX4_WS): skipping sando_px4_bridge"
    SKIP_ARGS=(--packages-skip sando_px4_bridge)
elif [ -f "$PX4_CACHE" ] && grep -q "^ament_cmake_core_DIR:PATH=/opt/ros/" "$PX4_CACHE" \
        && ! grep -q "^ament_cmake_core_DIR:PATH=/opt/ros/$ROS_DISTRO_ARG/" "$PX4_CACHE"; then
    # The distro px4_msgs was configured with (its cache can list several on the prefix path)
    other=$(grep "^ament_cmake_core_DIR:PATH=" "$PX4_CACHE" | cut -d/ -f4)
    echo "  px4_msgs in $PX4_WS was built for ROS 2 $other, not $ROS_DISTRO_ARG:"
    echo "  skipping sando_px4_bridge (rebuild px4_msgs with only $ROS_DISTRO_ARG sourced to enable it)"
    SKIP_ARGS=(--packages-skip sando_px4_bridge)
fi
# Livox driver (Gazebo Classic lidar simulation only), built by setup.sh
LIVOX_WS="$CODE_DIR/livox_ws"
if [ -d "$LIVOX_WS/install/livox_ros_driver2" ]; then
    export CMAKE_PREFIX_PATH="$LIVOX_WS/install/livox_ros_driver2:$CMAKE_PREFIX_PATH"
    export LD_LIBRARY_PATH="$LIVOX_WS/install/livox_ros_driver2/lib:$LD_LIBRARY_PATH"
fi

SELECT_ARGS=()
if [ -n "$PACKAGES" ]; then
    read -ra PKG_LIST <<< "$PACKAGES"
    SELECT_ARGS=(--packages-select "${PKG_LIST[@]}")
fi

# deps/gazebo_ros_pkgs overrides the /opt/ros copies on Humble; only colcon versions with the
# override check extension know --allow-overriding
OVERRIDE_ARGS=()
if colcon build --help 2>/dev/null | grep -q -- --allow-overriding; then
    OVERRIDE_ARGS=(--allow-overriding gazebo_dev gazebo_msgs gazebo_ros gazebo_ros_pkgs gazebo_plugins)
fi

cd "$SANDO_WS"
colcon build \
    --parallel-workers "$NUM_JOBS" \
    "${SELECT_ARGS[@]}" "${SKIP_ARGS[@]}" "${OVERRIDE_ARGS[@]}" \
    --cmake-args "${CMAKE_ARGS[@]}" -DSANDO_WITH_GUROBI="$WITH_GUROBI" -DSANDO_WITH_HIGHS_SHOT="$WITH_HIGHS_SHOT"
echo "$ROS_DISTRO_ARG" > "$SANDO_WS/install/.ros_distro"

echo ""
echo "============================================="
echo "  Build complete. To use it:"
echo "    source /opt/ros/$ROS_DISTRO_ARG/setup.bash"
echo "    source $DECOMP_WS/install/local_setup.bash"
echo "    source $SANDO_WS/install/local_setup.bash"
echo "============================================="
