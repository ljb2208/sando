FROM ros:humble

LABEL Maintainer="Kota Kondo <kkondo@mit.edu>" \
      Description="Ubuntu 22.04, ROS2 Humble for SANDO"

ENV ROS_DISTRO=humble
ENV DEBIAN_FRONTEND=noninteractive

# System packages
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -q -y --no-install-recommends \
        git tmux vim wget tmuxp make gdb && \
    rm -rf /var/lib/apt/lists/*

# ROS 2 dependencies
RUN apt-get update && apt-get install -y \
        ros-${ROS_DISTRO}-gazebo-ros-pkgs \
        ros-${ROS_DISTRO}-pcl-conversions \
        ros-${ROS_DISTRO}-example-interfaces \
        ros-${ROS_DISTRO}-pcl-ros \
        ros-${ROS_DISTRO}-rviz2 \
        ros-${ROS_DISTRO}-rviz-common \
        ros-${ROS_DISTRO}-rqt-gui \
        ros-${ROS_DISTRO}-rqt-gui-py \
        ros-${ROS_DISTRO}-tf2-tools \
        ros-${ROS_DISTRO}-tf-transformations \
        ros-${ROS_DISTRO}-desktop \
        ros-dev-tools \
        ros-${ROS_DISTRO}-gazebo-ros2-control \
        ros-${ROS_DISTRO}-xacro \
        libpcl-dev build-essential nlohmann-json3-dev && \
    rm -rf /var/lib/apt/lists/*

# Gurobi optimizer
RUN wget -q https://packages.gurobi.com/11.0/gurobi11.0.3_linux64.tar.gz -P /tmp && \
    tar -xzf /tmp/gurobi11.0.3_linux64.tar.gz -C /opt && \
    rm /tmp/gurobi11.0.3_linux64.tar.gz && \
    cd /opt/gurobi1103/linux64/src/build && make && cp libgurobi_c++.a ../../lib/
ENV GUROBI_HOME="/opt/gurobi1103/linux64"
ENV PATH="${PATH}:${GUROBI_HOME}/bin"
ENV LD_LIBRARY_PATH="${GUROBI_HOME}/lib"
# COPY gurobi.lic /opt/gurobi1103/

# Parallel build jobs, limited by default to avoid OOM kills (override with --build-arg BUILD_JOBS=N)
ARG BUILD_JOBS=2

# HiGHS (LP/QP/MIP solver) and SHOT (convex MINLP solver, with HiGHS as its MIP solver), built from
# their GitHub sources and installed to /usr/local. These back SANDO's default trajectory solver
# (SolverHighsShot). SHOT is built against this HiGHS so that SANDO and SHOT share one libhighs.
# Override the pinned versions with --build-arg HIGHS_VERSION=<tag> / --build-arg SHOT_COMMIT=<sha>.
ARG HIGHS_VERSION=v1.15.1
ARG SHOT_COMMIT=d9db115573f424e94b600ad59f704701e7ceded2
RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake pkg-config coinor-libipopt-dev && \
    rm -rf /var/lib/apt/lists/*
RUN git clone --branch ${HIGHS_VERSION} --depth 1 https://github.com/ERGO-Code/HiGHS.git /tmp/HiGHS && \
    cmake -S /tmp/HiGHS -B /tmp/HiGHS/build -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF && \
    cmake --build /tmp/HiGHS/build --parallel ${BUILD_JOBS} && \
    cmake --install /tmp/HiGHS/build && \
    rm -rf /tmp/HiGHS
RUN mkdir /tmp/SHOT && cd /tmp/SHOT && \
    git init -q && git remote add origin https://github.com/coin-or/SHOT.git && \
    git fetch -q --depth 1 origin ${SHOT_COMMIT} && git checkout -q FETCH_HEAD && \
    git submodule update --init --depth 1 ThirdParty/spdlog ThirdParty/CppAD ThirdParty/eigen && \
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local \
          -DHAS_HIGHS=ON -DUSE_EXTERNAL_HIGHS=ON -DHIGHS_EXTERNAL_DIR=/usr/local -DHAS_IPOPT=ON \
          -DHAS_CBC=OFF -DHAS_CPLEX=OFF -DHAS_GUROBI=OFF -DHAS_GAMS=OFF \
          -DBUILD_DOCUMENTATION=OFF && \
    cmake --build build --parallel ${BUILD_JOBS} && \
    cmake --install build && \
    ldconfig && \
    rm -rf /tmp/SHOT

# Clone SANDO from GitHub at a pinned tag (override with `make build SANDO_TAG=<tag>`).
# The default tracks `main`; tag a release and update SANDO_TAG here when publishing.
ENV SANDO_WS=/root/sando_ws
ENV LIVOX_WS=/root/livox_ws
RUN mkdir -p ${SANDO_WS}/src
WORKDIR ${SANDO_WS}/src
# DUMMY busts the cache of the clone below: `make build` sets it to the commit SHA of SANDO_TAG, so
# the clone (and everything after it) is redone exactly when new commits are pushed.
ARG DUMMY=unknown
#ARG SANDO_TAG=v0.0.3
ARG SANDO_TAG=main
ARG SANDO_REPO=https://github.com/ljb2208/sando.git
RUN DUMMY=${DUMMY} git clone --branch ${SANDO_TAG} --depth 1 ${SANDO_REPO} sando && \
    cd sando && \
    git submodule update --init --recursive

# Symlink deps for colcon package discovery (skip Livox-SDK2 and livox_ros_driver2 — built separately)
WORKDIR ${SANDO_WS}/src
RUN for d in sando/deps/*/; do \
      name=$(basename "$d"); \
      case "$name" in Livox-SDK2|livox_ros_driver2) continue;; esac; \
      [ ! -e "$name" ] && ln -s "$d" "$name"; \
    done

# Build Livox-SDK2 (cmake install to /usr/local)
WORKDIR ${SANDO_WS}/src/sando/deps/Livox-SDK2
RUN mkdir -p build && cd build && cmake .. && make -j && make install

# Build livox_ros_driver2 (in its own workspace — requires build.sh for package file setup)
RUN mkdir -p ${LIVOX_WS}/src
WORKDIR ${LIVOX_WS}/src
RUN cp -r ${SANDO_WS}/src/sando/deps/livox_ros_driver2 livox_ros_driver2 && \
    cd livox_ros_driver2 && \
    rm -f COLCON_IGNORE && \
    cp -f package_ROS2.xml package.xml && \
    cp -rf launch_ROS2/ launch/
SHELL ["/bin/bash", "-c"]
WORKDIR ${LIVOX_WS}
RUN source /opt/ros/humble/setup.bash && \
    colcon build --cmake-args -DCMAKE_BUILD_TYPE=Release \
                              -DROS_EDITION=ROS2 -DHUMBLE_ROS=humble

# Build all colcon packages (three passes: decomp_util, main, ros2_livox_simulation)
# Parallelism is limited to avoid OOM kills on template-heavy files like
# gazebo_ros_camera.cpp (which can use 3-4 GB RAM per compilation unit).
# Override with --build-arg BUILD_JOBS=N if you have more RAM available.
# SANDO uses HiGHS + SHOT by default; use --build-arg SANDO_USE_GUROBI=ON for Gurobi instead.
ARG SANDO_USE_GUROBI=OFF
ENV MAKEFLAGS="-j${BUILD_JOBS}"
WORKDIR ${SANDO_WS}
RUN source /opt/ros/humble/setup.bash && \
    colcon build --packages-select decomp_util \
                 --parallel-workers ${BUILD_JOBS} \
                 --cmake-args -DCMAKE_BUILD_TYPE=Release
RUN source /opt/ros/humble/setup.bash && \
    source ${SANDO_WS}/install/setup.bash && \
    colcon build --packages-ignore ros2_livox_simulation \
                 --parallel-workers ${BUILD_JOBS} \
                 --cmake-args -DCMAKE_BUILD_TYPE=Release -DSANDO_USE_GUROBI=${SANDO_USE_GUROBI}
RUN source /opt/ros/humble/setup.bash && \
    source ${SANDO_WS}/install/setup.bash && \
    source ${LIVOX_WS}/install/setup.bash && \
    colcon build --packages-select ros2_livox_simulation \
                 --parallel-workers ${BUILD_JOBS} \
                 --cmake-args -DCMAKE_BUILD_TYPE=Release

# Gazebo audio fix
RUN echo 'pcm.!default { type plug; slave.pcm "null" }' >> /usr/share/alsa/alsa.conf

# Livox library path
ENV LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}

# Shell aliases
RUN echo "alias tks='tmux kill-server'" >> ~/.bashrc && \
    echo "alias roscd='cd ${SANDO_WS}'" >> ~/.bashrc && \
    echo "alias ss='source ${SANDO_WS}/install/setup.bash'" >> ~/.bashrc && \
    echo "alias cbd='cd ${SANDO_WS} && colcon build --packages-select sando && ss'" >> ~/.bashrc && \
    echo "alias cbps='cd ${SANDO_WS} && colcon build --packages-select'" >> ~/.bashrc

# macOS visualization support via Xpra (browser-based RViz streaming)
# Used by `make run-mac-interactive` — dormant on Linux builds
RUN apt-get update && apt-get install -y --no-install-recommends \
    xpra xvfb libgl1-mesa-glx libgl1-mesa-dri psmisc \
    && rm -rf /var/lib/apt/lists/*

# Entry point
COPY sando.sh /
RUN chmod +x /sando.sh
ENTRYPOINT ["/sando.sh"]
