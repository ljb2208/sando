#!/usr/bin/env python3
"""Follow a list of GPS waypoints (geographic_msgs/GeoPath, e.g. from the AutoRally web OCS) by
handing SANDO one terminal goal at a time.

  waypoints     (geographic_msgs/GeoPath, latched)  lat/lon list; an empty list stops the vehicle
  local_origin  (sensor_msgs/NavSatFix, latched)    geodetic origin of the ENU frame the vehicle's
                                                    state is in (the AutoRally state estimator
                                                    publishes it), unless use_fixed_origin is set
  state         (dynus_interfaces/State)            vehicle position in that ENU frame
  sando_status  (std_msgs/String)                   SANDO's status (GOAL_REACHED at the last goal)
  ->
  term_goal       (geometry_msgs/PoseStamped)  the waypoint being driven to, in local x/y
  local_waypoints (nav_msgs/Path, latched)     the waypoints in local x/y, for RViz
  waypoint_status (std_msgs/String)            "idle", "waiting for ...", "following i/n", "done"

SANDO plans each goal to a stop, so the follower moves on to the next waypoint once the vehicle is
within switch_radius of the current one; only the last waypoint is driven to until SANDO reports
GOAL_REACHED (or the vehicle is within goal_radius of it).
"""

import math

import rclpy
from dynus_interfaces.msg import State
from geographic_msgs.msg import GeoPath
from geometry_msgs.msg import PoseStamped
from nav_msgs.msg import Path
from rclpy.node import Node
from rclpy.qos import DurabilityPolicy, QoSProfile, ReliabilityPolicy
from sensor_msgs.msg import NavSatFix
from std_msgs.msg import String

# WGS84
_A = 6378137.0
_F = 1.0 / 298.257223563
_E2 = _F * (2.0 - _F)


def geodetic_to_ecef(lat, lon, alt):
    phi, lam = math.radians(lat), math.radians(lon)
    sphi, cphi = math.sin(phi), math.cos(phi)
    n = _A / math.sqrt(1.0 - _E2 * sphi * sphi)
    return ((n + alt) * cphi * math.cos(lam),
            (n + alt) * cphi * math.sin(lam),
            (n * (1.0 - _E2) + alt) * sphi)


class LocalCartesian:
    """East/north/up about a geodetic origin: the same exact transformation as
    GeographicLib::LocalCartesian, which the AutoRally state estimator uses."""

    def __init__(self, lat0, lon0, alt0=0.0):
        self.origin = (lat0, lon0, alt0)
        self.x0, self.y0, self.z0 = geodetic_to_ecef(lat0, lon0, alt0)
        phi, lam = math.radians(lat0), math.radians(lon0)
        self.sphi, self.cphi = math.sin(phi), math.cos(phi)
        self.slam, self.clam = math.sin(lam), math.cos(lam)

    def forward(self, lat, lon, alt=None):
        if alt is None or not math.isfinite(alt):
            alt = self.origin[2]  # unknown altitude: on the origin's height
        x, y, z = geodetic_to_ecef(lat, lon, alt)
        dx, dy, dz = x - self.x0, y - self.y0, z - self.z0
        e = -self.slam * dx + self.clam * dy
        n = -self.sphi * self.clam * dx - self.sphi * self.slam * dy + self.cphi * dz
        u = self.cphi * self.clam * dx + self.cphi * self.slam * dy + self.sphi * dz
        return e, n, u


def thin(points, min_spacing):
    """Drop waypoints closer than min_spacing to the previous one kept. The last waypoint is always
    kept (it replaces the one before it when they are too close): the end of the route matters."""
    if not points:
        return [], 0
    out = [points[0]]
    for p in points[1:]:
        if math.dist(p, out[-1]) >= min_spacing:
            out.append(p)
    if out[-1] != points[-1]:
        out[-1] = points[-1]
    return out, len(points) - len(out)


def start_index(points, pos):
    """Waypoint to start from: the nearest one, or the one after it when the vehicle is already
    closer to that one than the nearest waypoint is (i.e. it is on the way there)."""
    i = min(range(len(points)), key=lambda k: math.dist(points[k], pos))
    if i + 1 < len(points) and math.dist(pos, points[i + 1]) < math.dist(points[i], points[i + 1]):
        i += 1
    return i


class WaypointFollower(Node):
    def __init__(self):
        super().__init__("waypoint_follower")
        p = lambda name, default: self.declare_parameter(name, default).value  # noqa: E731
        self.switch_radius = p("switch_radius", 2.0)  # [m] move on to the next waypoint
        self.goal_radius = p("goal_radius", 0.5)  # [m] last waypoint reached
        self.min_spacing = p("min_spacing", 1.0)  # [m] closer waypoints are dropped
        self.start_at_nearest = p("start_at_nearest", True)  # else always start at the first
        self.loop = p("loop", False)  # back to the first waypoint after the last
        self.goal_frame = p("goal_frame_id", "map")
        self.goal_z = p("default_goal_z", 0.0)  # SANDO's force_goal_z replaces it anyway
        use_fixed_origin = p("use_fixed_origin", False)
        origin = (p("origin_lat", 0.0), p("origin_lon", 0.0), p("origin_alt", 0.0))

        self.geo = []  # [(lat, lon, alt)] as received
        self.points = []  # [(x, y)] local, after thinning
        self.index = None  # waypoint being driven to
        self.done = False
        self.proj = None
        self.pos = None
        self.sando_status = ""

        latched = QoSProfile(depth=1, reliability=ReliabilityPolicy.RELIABLE,
                             durability=DurabilityPolicy.TRANSIENT_LOCAL)
        self.pub_goal = self.create_publisher(PoseStamped, "term_goal", 10)
        self.pub_path = self.create_publisher(Path, "local_waypoints", latched)
        self.pub_status = self.create_publisher(String, "waypoint_status", 10)

        self.create_subscription(GeoPath, p("waypoints_topic", "/waypoints"), self.on_waypoints, latched)
        if use_fixed_origin:
            self.set_origin(*origin, "parameters")
        else:
            self.create_subscription(NavSatFix, p("origin_topic", "/local_origin"), self.on_origin, latched)
        self.create_subscription(State, "state", self.on_state, 10)
        self.create_subscription(String, "sando_status", self.on_sando_status, 10)
        self.create_timer(0.1, self.step)

    # ------------------------------------------------------------------ inputs

    def set_origin(self, lat, lon, alt, source):
        if self.proj is not None and self.proj.origin == (lat, lon, alt):
            return
        self.proj = LocalCartesian(lat, lon, alt)
        self.get_logger().info(f"ENU origin ({source}): {lat:.9f}, {lon:.9f}, {alt:.2f} m")
        if self.geo:  # the local coordinates moved: convert again and pick up from the vehicle
            self.load(restart=True)

    def on_origin(self, msg):
        self.set_origin(msg.latitude, msg.longitude, msg.altitude, "local_origin")

    def on_waypoints(self, msg):
        self.geo = [(ps.pose.position.latitude, ps.pose.position.longitude, ps.pose.position.altitude)
                    for ps in msg.poses]
        self.get_logger().info(f"Received {len(self.geo)} waypoints")
        if not self.geo:
            self.points, self.index, self.done = [], None, False
            self.publish_path()
            self.stop()
            return
        self.load(restart=True)

    def on_state(self, msg):
        self.pos = (msg.pos.x, msg.pos.y)

    def on_sando_status(self, msg):
        self.sando_status = msg.data

    # ------------------------------------------------------------------ route

    def load(self, restart):
        if self.proj is None:
            return  # converted once the origin arrives
        local = [self.proj.forward(*g)[:2] for g in self.geo]
        self.points, dropped = thin(local, self.min_spacing)
        if dropped:
            self.get_logger().warn(
                f"Dropped {dropped} waypoint(s) closer than {self.min_spacing} m to the previous one")
        far = max(math.hypot(*q) for q in self.points)
        if far > 5000.0:
            self.get_logger().warn(f"A waypoint is {far / 1000:.1f} km from the ENU origin: wrong origin?")
        self.publish_path()
        if restart:
            self.index, self.done = None, False  # chosen in step() once the position is known

    def publish_path(self):
        path = Path()
        path.header.stamp = self.get_clock().now().to_msg()
        path.header.frame_id = self.goal_frame
        for (x, y) in self.points:
            ps = PoseStamped()
            ps.header = path.header
            ps.pose.position.x, ps.pose.position.y, ps.pose.position.z = x, y, self.goal_z
            ps.pose.orientation.w = 1.0
            path.poses.append(ps)
        self.pub_path.publish(path)

    def send_goal(self, i):
        self.index = i
        x, y = self.points[i]
        self.publish_goal(x, y)
        self.get_logger().info(f"Waypoint {i + 1}/{len(self.points)}: ({x:.1f}, {y:.1f})")

    def publish_goal(self, x, y):
        msg = PoseStamped()
        msg.header.stamp = self.get_clock().now().to_msg()
        msg.header.frame_id = self.goal_frame
        msg.pose.position.x, msg.pose.position.y, msg.pose.position.z = x, y, self.goal_z
        msg.pose.orientation.w = 1.0
        self.pub_goal.publish(msg)

    def stop(self):
        """SANDO has no cancel: a goal at the vehicle's position brings it to a stop there."""
        if self.pos is not None:
            self.publish_goal(*self.pos)
            self.get_logger().info("Waypoints cleared: stopping")

    def step(self):
        if not self.geo:
            self.set_status("idle")
        elif self.proj is None:
            self.set_status("waiting for local_origin")
        elif self.pos is None:
            self.set_status("waiting for state")
        elif self.pub_goal.get_subscription_count() == 0:
            self.set_status("waiting for SANDO (term_goal)")
        elif self.done:
            self.set_status("done")
        elif self.index is None:
            self.send_goal(start_index(self.points, self.pos) if self.start_at_nearest else 0)
        else:
            self.advance()
            if not self.done:
                self.set_status(f"following {self.index + 1}/{len(self.points)}")

    def advance(self):
        i = self.index
        # skip every waypoint already within reach, so close waypoints don't send a burst of goals
        while i + 1 < len(self.points) and math.dist(self.points[i], self.pos) < self.switch_radius:
            i += 1
        if i != self.index:
            self.send_goal(i)
            return
        if i + 1 < len(self.points):
            return
        d = math.dist(self.points[i], self.pos)
        if d < self.goal_radius or (self.sando_status == "GOAL_REACHED" and d < self.switch_radius):
            if self.loop and len(self.points) > 1:
                self.send_goal(0)
            else:
                self.done = True
                self.get_logger().info("Last waypoint reached")

    def set_status(self, text):
        self.pub_status.publish(String(data=text))


def main(args=None):
    rclpy.init(args=args)
    node = WaypointFollower()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.try_shutdown()


if __name__ == "__main__":
    main()
