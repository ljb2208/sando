#!/usr/bin/env python3

# /* ----------------------------------------------------------------------------
#  * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

"""Web interface for SANDO.

Subscribes to one agent's SANDO topics, reduces them (latest message only, voxelized clouds,
corridor wireframes), and serves a page that shows status, performance, parameters, obstacles and
the planned trajectory, updated live over Server-Sent Events. Parameter changes and goals are sent
back through ROS (the sando_node parameter services and term_goal).

Only the Python standard library, numpy and ROS are used, so it runs on the onboard computer as is:
    ros2 run sando sando_web.py --ros-args -p namespace:=NX01 -p port:=8080
and open http://<host>:8080 in a browser.
"""

import base64
import faulthandler
import itertools
import json
import os
import signal
import statistics
import threading
import time
from collections import deque
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse

import numpy as np
import rclpy
from rclpy.executors import MultiThreadedExecutor
from rclpy.node import Node
from rclpy.parameter import Parameter, parameter_value_to_python
from rclpy.parameter_client import AsyncParameterClient
from rclpy.qos import QoSProfile, ReliabilityPolicy, HistoryPolicy, qos_profile_sensor_data
from sensor_msgs.msg import PointCloud2
from sensor_msgs_py import point_cloud2
from geometry_msgs.msg import PointStamped, PoseStamped
from std_msgs.msg import Empty, String
from visualization_msgs.msg import Marker, MarkerArray
from decomp_ros_msgs.msg import PolyhedronArray
from dynus_interfaces.msg import ComputationTimes, State

try:
    from ament_index_python.packages import get_package_share_directory

    DEFAULT_WEB_ROOT = os.path.join(get_package_share_directory("sando"), "web")
except Exception:  # running from the source tree
    DEFAULT_WEB_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "web")

TIMING_FIELDS = [
    "total_replanning_ms",
    "global_planning_ms",
    "cvx_decomp_ms",
    "local_traj_ms",
    "safe_paths_ms",
    "safety_check_ms",
    "yaw_sequence_ms",
]

RELIABLE = QoSProfile(
    reliability=ReliabilityPolicy.RELIABLE, history=HistoryPolicy.KEEP_LAST, depth=10
)
# Visualization topics: only the latest message is shown
LATEST = QoSProfile(
    reliability=ReliabilityPolicy.RELIABLE, history=HistoryPolicy.KEEP_LAST, depth=1
)


def r3(v):
    return round(float(v), 3)


def color_hex(c):
    return (int(c.r * 255) << 16) | (int(c.g * 255) << 8) | int(c.b * 255)


def b64_f32(arr):
    return base64.b64encode(np.ascontiguousarray(arr, dtype=np.float32).tobytes()).decode()


# ---------------------------------------------------------------------------- geometry helpers


def markers_to_geometry(markers):
    """Flatten markers into line segments and points with per-element colors.

    Handles the marker types SANDO publishes (ARROW with two points for the colored trajectory,
    LINE_STRIP/LINE_LIST paths, SPHERE/CUBE lists and single shapes). Point lists are offset by
    the marker position; their orientation is ignored (SANDO publishes identity poses).
    """
    seg, seg_col, pts, pts_col = [], [], [], []
    size = 0.1
    for m in markers:
        if m.type == Marker.TEXT_VIEW_FACING or m.type == Marker.MESH_RESOURCE:
            continue
        ox, oy, oz = m.pose.position.x, m.pose.position.y, m.pose.position.z
        base = color_hex(m.color)

        def col(i):
            return color_hex(m.colors[i]) if i < len(m.colors) else base

        p = [(q.x + ox, q.y + oy, q.z + oz) for q in m.points]
        if m.type == Marker.ARROW and len(p) >= 2:
            seg += [*p[0], *p[1]]
            seg_col.append(base)
        elif m.type == Marker.LINE_STRIP:
            for i in range(len(p) - 1):
                seg += [*p[i], *p[i + 1]]
                seg_col.append(col(i))
        elif m.type == Marker.LINE_LIST:
            for i in range(0, len(p) - 1, 2):
                seg += [*p[i], *p[i + 1]]
                seg_col.append(col(i))
        elif m.type in (Marker.SPHERE_LIST, Marker.CUBE_LIST, Marker.POINTS):
            for i, q in enumerate(p):
                pts += q
                pts_col.append(col(i))
            size = max(size, m.scale.x)
        elif m.type in (Marker.SPHERE, Marker.CUBE, Marker.CYLINDER, Marker.ARROW):
            pts += [ox, oy, oz]
            pts_col.append(base)
            size = max(size, m.scale.x)
    return {
        "seg": b64_f32(seg) if seg else "",
        "segCol": seg_col,
        "pts": b64_f32(pts) if pts else "",
        "ptsCol": pts_col,
        "size": r3(size),
    }


def polyhedron_edges(points, normals):
    """Wireframe of the convex polyhedron {x : n_i . (x - p_i) <= 0}: returns a flat segment list.

    Vertices are the feasible intersections of plane triples; two vertices form an edge when they
    share at least two planes.
    """
    P = np.array([[p.x, p.y, p.z] for p in points])
    N = np.array([[n.x, n.y, n.z] for n in normals])
    k = len(N)
    if k < 4:
        return []
    b = np.einsum("ij,ij->i", N, P)
    tri = np.array(list(itertools.combinations(range(k), 3)))
    A = N[tri]  # (T, 3, 3)
    det = np.linalg.det(A)
    ok = np.abs(det) > 1e-9
    tri, A = tri[ok], A[ok]
    if len(tri) == 0:
        return []
    X = np.linalg.solve(A, b[tri][..., None])[..., 0]  # (T, 3)
    feasible = np.all(X @ N.T <= b + 1e-6, axis=1)
    X = X[feasible]
    if len(X) == 0:
        return []
    # Active planes of each (deduplicated) vertex
    keys = np.round(X, 4)
    _, idx = np.unique(keys, axis=0, return_index=True)
    V = X[idx]
    active = np.abs(V @ N.T - b) < 1e-5  # (V, k)
    shared = active.astype(np.int32) @ active.T.astype(np.int32)
    i, j = np.where(np.triu(shared >= 2, 1))
    return np.hstack([V[i], V[j]]).ravel().tolist()


def voxelize(xyz, voxel, max_points):
    if len(xyz) == 0:
        return xyz
    if voxel > 0:
        q = np.floor(xyz / voxel).astype(np.int64)
        _, idx = np.unique(q, axis=0, return_index=True)
        xyz = (q[idx] + 0.5) * voxel
    if len(xyz) > max_points:
        xyz = xyz[np.random.default_rng(0).choice(len(xyz), max_points, replace=False)]
    return xyz


# ---------------------------------------------------------------------------- ROS side


class SandoWeb(Node):
    def __init__(self):
        super().__init__("sando_web")
        self.ns = self.declare_parameter("namespace", "NX01").value.strip("/")
        self.port = self.declare_parameter("port", 8080).value
        self.host = self.declare_parameter("host", "0.0.0.0").value
        self.web_root = self.declare_parameter("web_root", DEFAULT_WEB_ROOT).value
        self.read_only = self.declare_parameter("read_only", False).value
        self.voxel = self.declare_parameter("cloud_voxel_size", 0.2).value
        self.max_cloud_points = self.declare_parameter("max_cloud_points", 30000).value
        self.cloud_topics = self.declare_parameter(
            "cloud_topics",
            ["occupancy_grid", "sando_occupied_cloud", "/map_generator/global_cloud"],
        ).value
        self.marker_topics = self.declare_parameter(
            "marker_topics",
            ["traj_committed_colored", "hgp_path_marker", "cp", "actual_traj"],
        ).value
        self.goal_z = self.declare_parameter("default_goal_z", 3.0).value

        ns = f"/{self.ns}"
        self.sando_node_name = f"{ns}/sando_node"
        self.lock = threading.Lock()

        # Latest raw data (converted in the broadcaster thread, not in the callbacks)
        self.state = None
        self.trail = deque(maxlen=1500)
        self.last_trail_t = 0.0
        self.status = "UNKNOWN"
        self.status_t = 0.0
        self.points = {}  # A, G, E, G_term -> [x, y, z]
        self.term_goal = None
        self.markers = {t: {} for t in self.marker_topics}  # topic -> {(ns, id): Marker}
        self.markers_dirty = {t: False for t in self.marker_topics}
        self.marker_geom = {t: None for t in self.marker_topics}
        self.poly_msg = None
        self.poly_edges = []
        self.poly_dirty = False
        self.clouds = {t: None for t in self.cloud_topics}
        self.clouds_dirty = {t: False for t in self.cloud_topics}
        self.obstacles = []  # latest tracked_obstacles boxes
        self.obstacles_t = 0.0
        self.timings = deque(maxlen=8000)  # ~60 s at the 100 Hz replan rate  # (wall time, success, {field: ms}, factor)
        self.events = deque(maxlen=200)
        self.goals_reached = 0

        self.create_subscription(State, f"{ns}/state", self.on_state, qos_profile_sensor_data)
        self.create_subscription(String, f"{ns}/sando_status", self.on_status, RELIABLE)
        self.create_subscription(
            ComputationTimes, f"{ns}/computation_times", self.on_times, QoSProfile(depth=100)
        )
        self.create_subscription(Empty, f"{ns}/goal_reached", self.on_goal_reached, RELIABLE)
        self.create_subscription(PoseStamped, f"{ns}/term_goal", self.on_term_goal, RELIABLE)
        for name in ("A", "G", "E", "G_term"):
            self.create_subscription(
                PointStamped,
                f"{ns}/point_{name}",
                lambda m, n=name: self.on_point(n, m),
                LATEST,
            )
        for t in self.marker_topics:
            self.create_subscription(
                MarkerArray, self.resolve(t), lambda m, t=t: self.on_markers(t, m), LATEST
            )
        self.create_subscription(PolyhedronArray, f"{ns}/poly_safe", self.on_poly, LATEST)
        for t in self.cloud_topics:
            self.create_subscription(
                PointCloud2, self.resolve(t), lambda m, t=t: self.on_cloud(t, m),
                qos_profile_sensor_data,
            )
        # The obstacles SANDO tracks, as SANDO publishes them at 5 Hz. Not /trajs itself: it carries
        # one message per obstacle update (> 1 kHz in the dynamic sims), too many for Python
        self.create_subscription(
            MarkerArray, f"{ns}/tracked_obstacles", self.on_obstacles, LATEST
        )

        self.pub_goal = self.create_publisher(PoseStamped, f"{ns}/term_goal", RELIABLE)
        self.param_client = AsyncParameterClient(self, self.sando_node_name)

    def resolve(self, topic):
        return topic if topic.startswith("/") else f"/{self.ns}/{topic}"

    def log_event(self, kind, text):
        with self.lock:
            self.events.append({"t": time.time(), "kind": kind, "text": text})

    # ------------------------------------------------------------------ callbacks

    def on_state(self, m):
        now = time.time()
        with self.lock:
            self.state = m
            if now - self.last_trail_t > 0.1:
                self.trail.append((m.pos.x, m.pos.y, m.pos.z))
                self.last_trail_t = now

    def on_status(self, m):
        changed = False
        with self.lock:
            if m.data != self.status:
                changed = self.status != "UNKNOWN" or m.data != "GOAL_REACHED"
                prev, self.status = self.status, m.data
            self.status_t = time.time()
        if changed:
            self.log_event("status", f"{prev} → {m.data}")

    def on_times(self, m):
        # SANDO publishes this only for replans it attempted (not while idle at the goal, ...)
        with self.lock:
            self.timings.append(
                (time.time(), m.result, {f: getattr(m, f) for f in TIMING_FIELDS},
                 m.successful_factor)
            )

    def on_goal_reached(self, _):
        with self.lock:
            self.goals_reached += 1
        self.log_event("goal", "Goal reached")

    def on_term_goal(self, m):
        p = m.pose.position
        with self.lock:
            self.term_goal = [r3(p.x), r3(p.y), r3(p.z)]
        self.log_event("goal", f"New goal ({p.x:.1f}, {p.y:.1f}, {p.z:.1f})")

    def on_point(self, name, m):
        with self.lock:
            self.points[name] = [r3(m.point.x), r3(m.point.y), r3(m.point.z)]

    def on_markers(self, topic, msg):
        with self.lock:
            store = self.markers[topic]
            for mk in msg.markers:
                if mk.action == Marker.DELETEALL:
                    store.clear()
                elif mk.action == Marker.DELETE:
                    store.pop((mk.ns, mk.id), None)
                else:
                    store[(mk.ns, mk.id)] = mk
            self.markers_dirty[topic] = True

    def on_poly(self, msg):
        with self.lock:
            self.poly_msg = msg
            self.poly_dirty = True

    def on_cloud(self, topic, msg):
        with self.lock:
            self.clouds[topic] = msg
            self.clouds_dirty[topic] = True

    def on_obstacles(self, msg):
        boxes = []
        for mk in msg.markers:
            if mk.action != Marker.ADD:
                continue
            p, s = mk.pose.position, mk.scale
            boxes.append({"id": mk.id, "agent": mk.ns == "agents",
                          "pos": [r3(p.x), r3(p.y), r3(p.z)], "bbox": [r3(s.x), r3(s.y), r3(s.z)]})
        with self.lock:
            self.obstacles = boxes
            self.obstacles_t = time.time()

    # ------------------------------------------------------------------ snapshots for the page

    def snap_state(self):
        with self.lock:
            s = self.state
            status, status_age = self.status, time.time() - self.status_t
            goal, pts = self.term_goal, dict(self.points)
        if s is None:
            return {"connected": status_age < 2.0, "status": status}
        pos = [r3(s.pos.x), r3(s.pos.y), r3(s.pos.z)]
        vel = [r3(s.vel.x), r3(s.vel.y), r3(s.vel.z)]
        q = s.quat
        out = {
            "connected": status_age < 2.0,
            "status": status,
            "pos": pos,
            "vel": vel,
            "quat": [r3(q.x), r3(q.y), r3(q.z), r3(q.w)],
            "speed": r3(np.linalg.norm(vel)),
            "goal": goal,
            "points": pts,
        }
        if goal:
            out["goalDist"] = r3(np.linalg.norm(np.subtract(goal, pos)))
        return out

    def snap_scene(self):
        now = time.time()
        with self.lock:
            dirty = [t for t, d in self.markers_dirty.items() if d]
            stores = {t: list(self.markers[t].values()) for t in dirty}
            for t in dirty:
                self.markers_dirty[t] = False
            trail = list(self.trail)
            obstacles = self.obstacles if now - self.obstacles_t < 2.0 else []
        for t, mks in stores.items():
            self.marker_geom[t] = markers_to_geometry(mks)
        return {
            "markers": {t: g for t, g in self.marker_geom.items() if g is not None},
            "trail": b64_f32(np.array(trail).ravel()) if trail else "",
            "obstacles": obstacles,
        }

    def snap_polys(self):
        with self.lock:
            msg, dirty = self.poly_msg, self.poly_dirty
            self.poly_dirty = False
        if msg is None:
            return None
        if dirty:
            edges = []
            for poly in msg.polyhedrons:
                try:
                    edges += polyhedron_edges(poly.points, poly.normals)
                except np.linalg.LinAlgError:
                    pass
            self.poly_edges = edges
        return {"seg": b64_f32(self.poly_edges) if self.poly_edges else "",
                "count": len(msg.polyhedrons)}

    def snap_clouds(self):
        """Only the clouds that changed since the last call."""
        with self.lock:
            todo = {t: self.clouds[t] for t, d in self.clouds_dirty.items() if d}
            for t in todo:
                self.clouds_dirty[t] = False
        out = {}
        for t, msg in todo.items():
            try:
                xyz = point_cloud2.read_points_numpy(msg, field_names=("x", "y", "z"),
                                                     skip_nans=True)
            except Exception as e:  # noqa: BLE001 - unexpected cloud layout
                self.get_logger().warn(f"cannot read cloud {t}: {e}", throttle_duration_sec=10.0)
                continue
            xyz = voxelize(np.asarray(xyz, dtype=np.float64), self.voxel, self.max_cloud_points)
            out[t] = {"pts": b64_f32(xyz.ravel()), "n": len(xyz), "voxel": self.voxel,
                      "frame": msg.header.frame_id}
        return out

    def snap_perf(self, window_s=60.0):
        now = time.time()
        with self.lock:
            rows = [r for r in self.timings if now - r[0] <= window_s]
            reached = self.goals_reached
        recent = [r for r in rows if now - r[0] <= 5.0]
        last100 = rows[-100:]
        ok = [r for r in last100 if r[1]]

        def stats(field):
            vals = [r[2][field] for r in ok if r[2][field] > 0]
            if not vals:
                return None
            vals.sort()
            return {
                "median": r3(statistics.median(vals)),
                "p95": r3(vals[min(len(vals) - 1, int(0.95 * len(vals)))]),
                "max": r3(vals[-1]),
            }

        series = {
            "t": [r3(r[0] - now) for r in rows],
            "ok": [1 if r[1] else 0 for r in rows],
            "total": [r3(r[2]["total_replanning_ms"]) for r in rows],
            "local": [r3(r[2]["local_traj_ms"]) for r in rows],
            "global": [r3(r[2]["global_planning_ms"]) for r in rows],
            "decomp": [r3(r[2]["cvx_decomp_ms"]) for r in rows],
        }
        # Decimate the chart series to keep updates small
        step = max(1, len(rows) // 600)
        series = {k: v[::step] for k, v in series.items()}
        return {
            "replanHz": r3(len(recent) / 5.0),
            "successPct": r3(100.0 * len(ok) / len(last100)) if last100 else None,
            "successHz": r3(sum(1 for r in recent if r[1]) / 5.0),
            "samples": len(last100),
            "goalsReached": reached,
            "lastFactor": r3(ok[-1][3]) if ok else None,
            "stats": {f: stats(f) for f in TIMING_FIELDS},
            "series": series,
        }

    def snap_events(self, since):
        with self.lock:
            return [e for e in self.events if e["t"] > since]

    # ------------------------------------------------------------------ parameters / goal

    def _wait(self, future, timeout=3.0):
        done = threading.Event()
        future.add_done_callback(lambda _: done.set())
        if not done.wait(timeout):
            raise TimeoutError(f"{self.sando_node_name} did not answer (is SANDO running?)")
        return future.result()

    def get_params(self):
        if not self.param_client.wait_for_services(timeout_sec=1.0):
            raise TimeoutError(f"no parameter services on {self.sando_node_name}")
        names = sorted(self._wait(self.param_client.list_parameters()).result.names)
        values = self._wait(self.param_client.get_parameters(names)).values
        params = {n: parameter_value_to_python(v) for n, v in zip(names, values)}
        live = set(params.pop("live_parameters", []) or [])
        out = []
        for n, v in params.items():
            if n == "use_sim_time" or n.startswith("qos_overrides"):
                continue
            if isinstance(v, (list, tuple)) or hasattr(v, "tolist"):
                v = list(v.tolist() if hasattr(v, "tolist") else v)
            out.append({"name": n, "value": v, "live": n in live,
                        "type": type(v).__name__ if not isinstance(v, list) else "list"})
        return {"node": self.sando_node_name, "params": out, "readOnly": self.read_only}

    def set_param(self, name, value):
        if self.read_only:
            return False, "the web interface is read-only (read_only:=true)"
        current = self._wait(self.param_client.get_parameters([name])).values[0]
        cur = parameter_value_to_python(current)
        try:
            if isinstance(cur, bool):
                v = value if isinstance(value, bool) else str(value).lower() in ("1", "true", "yes")
            elif isinstance(cur, int):
                v = int(value)
            elif isinstance(cur, float):
                v = float(value)
            elif isinstance(cur, str):
                v = str(value)
            else:
                return False, f"'{name}' has a type the web interface cannot set"
        except (TypeError, ValueError):
            return False, f"invalid value for '{name}': {value!r}"
        res = self._wait(self.param_client.set_parameters([Parameter(name, value=v)]))
        r = res.results[0]
        if r.successful:
            self.log_event("param", f"{name} = {v}")
        return r.successful, r.reason

    def send_goal(self, x, y, z):
        msg = PoseStamped()
        msg.header.frame_id = "map"
        msg.header.stamp = self.get_clock().now().to_msg()
        msg.pose.position.x, msg.pose.position.y = float(x), float(y)
        msg.pose.position.z = float(z if z is not None else self.goal_z)
        msg.pose.orientation.w = 1.0
        self.pub_goal.publish(msg)


# ---------------------------------------------------------------------------- web side


class Client:
    """Updates waiting to be written to one page: the latest of each kind only (event lists are
    kept in full), so a slow page skips intermediate updates instead of losing whole kinds."""

    ORDER = ("clouds", "polys", "scene", "perf", "state", "events")

    def __init__(self):
        self.cv = threading.Condition()
        self.pending = {}  # key -> payload; key is the kind, or "clouds:<topic>"
        self.events = []

    def put(self, key, payload):
        with self.cv:
            if key == "events":
                self.events.append(payload)
            else:
                self.pending[key] = payload
            self.cv.notify()

    def take(self, timeout):
        """Payloads to write, in a fixed order (empty after `timeout` seconds without updates)."""
        with self.cv:
            self.cv.wait_for(lambda: self.pending or self.events, timeout)
            pending, events = self.pending, self.events
            self.pending, self.events = {}, []
        rank = {k: i for i, k in enumerate(self.ORDER)}
        keys = sorted(pending, key=lambda k: rank.get(k.split(":")[0], len(rank)))
        return [pending[k] for k in keys] + events


class Broadcaster:
    """Builds the updates at fixed rates and fans them out to every connected page."""

    RATES = {"state": 15.0, "scene": 5.0, "polys": 2.0, "perf": 2.0, "clouds": 1.0, "events": 2.0}

    def __init__(self, node):
        self.node = node
        self.clients = set()
        self.lock = threading.Lock()
        self.last_event_t = 0.0
        self.full_clouds = {}  # latest encoded cloud per topic, sent to new pages
        self.full_polys = None
        threading.Thread(target=self.run, daemon=True).start()

    def add(self):
        c = Client()
        with self.lock:
            self.clients.add(c)
            for t, cloud in self.full_clouds.items():
                c.put(f"clouds:{t}", self.encode("clouds", {t: cloud}))
            if self.full_polys:
                c.put("polys", self.encode("polys", self.full_polys))
        return c

    def remove(self, c):
        with self.lock:
            self.clients.discard(c)

    @staticmethod
    def encode(kind, data):
        return f"event: {kind}\ndata: {json.dumps(data, separators=(',', ':'))}\n\n".encode()

    def send(self, kind, data):
        if kind == "clouds":  # one entry per topic, so a slow page still gets every cloud
            for t, cloud in data.items():
                self.send_payload(f"clouds:{t}", self.encode("clouds", {t: cloud}))
        else:
            self.send_payload(kind, self.encode(kind, data))

    def send_payload(self, key, payload):
        with self.lock:
            for c in list(self.clients):
                c.put(key, payload)

    def run(self):
        next_t = {k: 0.0 for k in self.RATES}
        makers = {
            "state": self.node.snap_state,
            "scene": self.node.snap_scene,
            "polys": self.node.snap_polys,
            "perf": self.node.snap_perf,
            "clouds": self.node.snap_clouds,
        }
        while rclpy.ok():
            now = time.time()
            with self.lock:
                idle = not self.clients
            for kind, rate in self.RATES.items():
                if now < next_t[kind]:
                    continue
                next_t[kind] = now + 1.0 / rate
                try:
                    if kind == "events":
                        ev = self.node.snap_events(self.last_event_t)
                        if ev:
                            self.last_event_t = ev[-1]["t"]
                            self.send("events", ev)
                        continue
                    if idle and kind != "clouds" and kind != "polys":
                        continue
                    data = makers[kind]()
                    if kind == "clouds":
                        self.full_clouds.update(data)
                        if not data:
                            continue
                    if kind == "polys":
                        if data is None:
                            continue
                        self.full_polys = data
                    self.send(kind, data)
                except Exception as e:  # noqa: BLE001 - keep serving
                    self.node.get_logger().error(f"{kind} update failed: {e}",
                                                 throttle_duration_sec=5.0)
            time.sleep(0.01)


def make_handler(node, broadcaster):
    root = os.path.realpath(node.web_root)
    types = {".html": "text/html", ".js": "text/javascript", ".css": "text/css",
             ".svg": "image/svg+xml", ".txt": "text/plain"}

    class Handler(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"

        def log_message(self, fmt, *args):  # quiet
            pass

        def reply(self, code, body, ctype="application/json"):
            if not isinstance(body, bytes):
                body = json.dumps(body).encode()
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            path = urlparse(self.path).path
            if path == "/events":
                return self.stream()
            if path == "/api/params":
                try:
                    return self.reply(200, node.get_params())
                except Exception as e:  # noqa: BLE001
                    return self.reply(503, {"error": str(e)})
            if path == "/api/info":
                return self.reply(200, {
                    "namespace": node.ns, "node": node.sando_node_name,
                    "clouds": node.cloud_topics, "markers": node.marker_topics,
                    "readOnly": node.read_only, "defaultGoalZ": node.goal_z,
                })
            rel = "index.html" if path in ("/", "") else path.lstrip("/")
            full = os.path.realpath(os.path.join(root, rel))
            if not full.startswith(root + os.sep) or not os.path.isfile(full):
                return self.reply(404, b"not found", "text/plain")
            with open(full, "rb") as f:
                self.reply(200, f.read(), types.get(os.path.splitext(full)[1], "application/octet-stream"))

        def do_POST(self):
            path = urlparse(self.path).path
            try:
                body = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
            except ValueError:
                return self.reply(400, {"ok": False, "reason": "invalid JSON"})
            try:
                if path == "/api/params":
                    ok, reason = node.set_param(body["name"], body["value"])
                    return self.reply(200, {"ok": ok, "reason": reason})
                if path == "/api/goal":
                    if node.read_only:
                        return self.reply(403, {"ok": False, "reason": "read-only"})
                    node.send_goal(body["x"], body["y"], body.get("z"))
                    return self.reply(200, {"ok": True})
            except KeyError as e:
                return self.reply(400, {"ok": False, "reason": f"missing {e}"})
            except Exception as e:  # noqa: BLE001
                return self.reply(503, {"ok": False, "reason": str(e)})
            self.reply(404, {"ok": False, "reason": "not found"})

        def stream(self):
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Connection", "keep-alive")
            self.end_headers()
            client = broadcaster.add()
            try:
                self.wfile.write(b"retry: 2000\n\n")
                while rclpy.ok():
                    payloads = client.take(timeout=10.0) or [b": keepalive\n\n"]
                    self.wfile.write(b"".join(payloads))
                    self.wfile.flush()
            except (BrokenPipeError, ConnectionResetError, OSError):
                pass
            finally:
                broadcaster.remove(client)
                self.close_connection = True

    return Handler


def main(args=None):
    # kill -USR1 <pid> prints every thread's stack (diagnosing a stuck update loop)
    faulthandler.register(signal.SIGUSR1, all_threads=True)
    rclpy.init(args=args)
    node = SandoWeb()
    # The events executor costs much less per callback than the multithreaded one
    try:
        from rclpy.experimental import EventsExecutor

        executor = EventsExecutor()
    except ImportError:
        executor = MultiThreadedExecutor(num_threads=4)
    executor.add_node(node)
    threading.Thread(target=executor.spin, daemon=True).start()

    broadcaster = Broadcaster(node)
    server = ThreadingHTTPServer((node.host, node.port), make_handler(node, broadcaster))
    server.daemon_threads = True
    shown = "localhost" if node.host in ("0.0.0.0", "") else node.host
    node.get_logger().info(
        f"SANDO web interface for /{node.ns}: http://{shown}:{node.port} (serving {node.web_root})"
    )
    print(f"SANDO web interface: http://{shown}:{node.port}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        executor.shutdown()
        node.destroy_node()
        if rclpy.ok():
            rclpy.shutdown()


if __name__ == "__main__":
    main()
