// SANDO web interface: live view of one agent's planner over Server-Sent Events (/events).
import * as THREE from "three";
import { OrbitControls } from "./vendor/OrbitControls.js";

const $ = (id) => document.getElementById(id);
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const fmt = (v, d = 1) => (v === null || v === undefined || Number.isNaN(v) ? "—" : Number(v).toFixed(d));

function f32(b64) {
  if (!b64) return new Float32Array(0);
  const bin = atob(b64);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return new Float32Array(bytes.buffer);
}

let info = { namespace: "", clouds: [], markers: [], readOnly: false, defaultGoalZ: 3.0 };

// ============================================================================ 3D scene

const container = $("scene");
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
container.appendChild(renderer.domElement);

const scene = new THREE.Scene();
THREE.Object3D.DEFAULT_UP.set(0, 0, 1);
const camera = new THREE.PerspectiveCamera(55, 1, 0.05, 2000);
camera.up.set(0, 0, 1);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
function resetView() {
  const p = vehicle.position;
  camera.position.set(p.x - 12, p.y - 12, p.z + 10);
  controls.target.copy(p);
}

scene.add(new THREE.AmbientLight(0xffffff, 0.8));
const sun = new THREE.DirectionalLight(0xffffff, 0.8);
sun.position.set(10, -10, 20);
scene.add(sun);

let grid = null;
function applyTheme() {
  scene.background = new THREE.Color(css("--scene-bg"));
  if (grid) scene.remove(grid);
  grid = new THREE.GridHelper(200, 200, css("--axis"), css("--scene-grid"));
  grid.rotation.x = Math.PI / 2; // XY plane, z up
  grid.material.transparent = true;
  grid.material.opacity = 0.6;
  scene.add(grid);
}
applyTheme();
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => {
  applyTheme();
  chart.draw();
});

// Layers: name -> { group, label, swatch, on }
const layers = new Map();
function layer(key, label, swatchColor, isPoint = false, on = true) {
  if (layers.has(key)) return layers.get(key).group;
  const group = new THREE.Group();
  group.visible = on;
  scene.add(group);
  layers.set(key, { group, label, swatchColor, isPoint });
  renderLayerList();
  return group;
}
function renderLayerList() {
  const list = $("layer-list");
  list.innerHTML = "";
  for (const [key, l] of layers) {
    const row = document.createElement("label");
    const cb = document.createElement("input");
    cb.type = "checkbox";
    cb.checked = l.group.visible;
    cb.addEventListener("change", () => { l.group.visible = cb.checked; });
    const sw = document.createElement("span");
    sw.className = "swatch" + (l.isPoint ? " pt" : "");
    sw.style.background = l.swatchColor;
    row.append(cb, sw, document.createTextNode(l.label));
    row.dataset.key = key;
    list.appendChild(row);
  }
}
function clearGroup(g) {
  for (const c of [...g.children]) {
    g.remove(c);
    c.geometry?.dispose();
    c.material?.dispose();
  }
}

// Vehicle: body + heading arrow, oriented by the state quaternion
const vehicle = new THREE.Group();
{
  const body = new THREE.Mesh(new THREE.BoxGeometry(0.5, 0.5, 0.12),
    new THREE.MeshStandardMaterial({ color: 0x3987e5 }));
  vehicle.add(body);
  vehicle.add(new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), new THREE.Vector3(), 0.8, 0xeb6834, 0.2, 0.12));
  vehicle.add(new THREE.AxesHelper(0.4));
}
layer("vehicle", "Vehicle", "#3987e5", true).add(vehicle);

const trailLine = new THREE.Line(new THREE.BufferGeometry(),
  new THREE.LineBasicMaterial({ color: 0x9085e9 }));
layer("trail", "Flown path", "#9085e9").add(trailLine);

const goalGroup = layer("goal", "Goal (G term), A / G / E", "#e34948", true);
const pointColors = { goal: 0xe34948, A: 0x1baf7a, G: 0xeda100, E: 0xe87ba4 };
const pointMeshes = {};
for (const [k, c] of Object.entries(pointColors)) {
  const m = new THREE.Mesh(new THREE.SphereGeometry(k === "goal" ? 0.3 : 0.15, 16, 12),
    new THREE.MeshStandardMaterial({ color: c }));
  m.visible = false;
  goalGroup.add(m);
  pointMeshes[k] = m;
}

const markerLabels = {
  traj_committed_colored: ["Planned trajectory", "#1baf7a"],
  hgp_path_marker: ["Global path", "#eda100"],
  cp: ["Control points", "#e87ba4", true],
  actual_traj: ["Actual trajectory", "#6250d6"],
};
function updateMarkers(markers) {
  for (const [topic, g] of Object.entries(markers)) {
    const [label, sw, pt] = markerLabels[topic] || [topic, "#898781"];
    const group = layer(`m:${topic}`, label, sw, !!pt, topic !== "actual_traj");
    clearGroup(group);
    const seg = f32(g.seg);
    if (seg.length) {
      const colors = new Float32Array(seg.length);
      const c = new THREE.Color();
      for (let i = 0; i < g.segCol.length; i++) {
        c.setHex(g.segCol[i]);
        for (let v = 0; v < 2; v++) colors.set([c.r, c.g, c.b], i * 6 + v * 3);
      }
      const geo = new THREE.BufferGeometry();
      geo.setAttribute("position", new THREE.BufferAttribute(seg, 3));
      geo.setAttribute("color", new THREE.BufferAttribute(colors, 3));
      group.add(new THREE.LineSegments(geo, new THREE.LineBasicMaterial({ vertexColors: true })));
    }
    const pts = f32(g.pts);
    if (pts.length) {
      const colors = new Float32Array(pts.length);
      const c = new THREE.Color();
      g.ptsCol.forEach((h, i) => { c.setHex(h); colors.set([c.r, c.g, c.b], i * 3); });
      const geo = new THREE.BufferGeometry();
      geo.setAttribute("position", new THREE.BufferAttribute(pts, 3));
      geo.setAttribute("color", new THREE.BufferAttribute(colors, 3));
      group.add(new THREE.Points(geo, new THREE.PointsMaterial({ size: Math.max(g.size, 0.08), vertexColors: true })));
    }
  }
}

const polyGroup = layer("polys", "Safe corridors", "#199e70", false, true);
function updatePolys(p) {
  clearGroup(polyGroup);
  const seg = f32(p.seg);
  if (!seg.length) return;
  const geo = new THREE.BufferGeometry();
  geo.setAttribute("position", new THREE.BufferAttribute(seg, 3));
  polyGroup.add(new THREE.LineSegments(geo,
    new THREE.LineBasicMaterial({ color: 0x1baf7a, transparent: true, opacity: 0.45 })));
}

// Clouds: one layer per topic, colored by height (sequential blue ramp)
const ramp = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"].map((h) => new THREE.Color(h));
function heightColor(t, out) {
  const x = Math.min(Math.max(t, 0), 1) * (ramp.length - 1);
  const i = Math.min(Math.floor(x), ramp.length - 2);
  return out.copy(ramp[i]).lerp(ramp[i + 1], x - i);
}
const cloudLabels = {
  occupancy_grid: "Occupied map (mapper)",
  sando_occupied_cloud: "Occupied map (SANDO)",
  "/map_generator/global_cloud": "Static obstacles (sim)",
};
function updateClouds(clouds) {
  for (const [topic, c] of Object.entries(clouds)) {
    const group = layer(`c:${topic}`, `${cloudLabels[topic] || topic} (${c.n})`, "#3987e5", true);
    layers.get(`c:${topic}`).label = `${cloudLabels[topic] || topic} (${c.n})`;
    clearGroup(group);
    const pts = f32(c.pts);
    if (pts.length) {
      let zmin = Infinity, zmax = -Infinity;
      for (let i = 2; i < pts.length; i += 3) { zmin = Math.min(zmin, pts[i]); zmax = Math.max(zmax, pts[i]); }
      const colors = new Float32Array(pts.length);
      const col = new THREE.Color();
      for (let i = 0; i < pts.length; i += 3) {
        heightColor((pts[i + 2] - zmin) / Math.max(zmax - zmin, 1e-3), col);
        colors[i] = col.r; colors[i + 1] = col.g; colors[i + 2] = col.b;
      }
      const geo = new THREE.BufferGeometry();
      geo.setAttribute("position", new THREE.BufferAttribute(pts, 3));
      geo.setAttribute("color", new THREE.BufferAttribute(colors, 3));
      group.add(new THREE.Points(geo, new THREE.PointsMaterial({ size: c.voxel || 0.15, vertexColors: true })));
    }
    renderLayerList();
  }
}

const obstGroup = layer("obst", "Dynamic obstacles / other agents", "#d03b3b");
const obstMat = new THREE.MeshStandardMaterial({ color: 0xd03b3b, transparent: true, opacity: 0.35 });
const agentMat = new THREE.MeshStandardMaterial({ color: 0x9085e9, transparent: true, opacity: 0.5 });
function updateObstacles(list) {
  for (const c of [...obstGroup.children]) { obstGroup.remove(c); c.geometry.dispose(); }
  for (const o of list) {
    const m = new THREE.Mesh(new THREE.BoxGeometry(...o.bbox), o.agent ? agentMat : obstMat);
    m.position.set(...o.pos);
    obstGroup.add(m);
  }
}

function updateTrail(b64) {
  const pts = f32(b64);
  trailLine.geometry.dispose();
  trailLine.geometry = new THREE.BufferGeometry();
  trailLine.geometry.setAttribute("position", new THREE.BufferAttribute(pts, 3));
}

let followVehicle = true;
let viewInitialized = false;
$("follow").addEventListener("change", (e) => { followVehicle = e.target.checked; });
$("reset-view").addEventListener("click", resetView);
if (window.matchMedia("(max-width: 900px)").matches) $("layers").open = false;

function resize() {
  const w = container.clientWidth, h = container.clientHeight;
  renderer.setSize(w, h);
  camera.aspect = w / Math.max(h, 1);
  camera.updateProjectionMatrix();
}
new ResizeObserver(resize).observe(container);
resize();
resetView();

const lastVehiclePos = new THREE.Vector3();
(function animate() {
  requestAnimationFrame(animate);
  if (followVehicle) {
    const delta = vehicle.position.clone().sub(lastVehiclePos);
    camera.position.add(delta);
    controls.target.add(delta);
  }
  lastVehiclePos.copy(vehicle.position);
  controls.update();
  renderer.render(scene, camera);
})();

// Shift+click the ground (z = 0 plane) to send a goal there at the default goal height
renderer.domElement.addEventListener("click", (e) => {
  if (!e.shiftKey) return;
  const r = renderer.domElement.getBoundingClientRect();
  const ndc = new THREE.Vector2(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
  const ray = new THREE.Raycaster();
  ray.setFromCamera(ndc, camera);
  const hit = new THREE.Vector3();
  if (ray.ray.intersectPlane(new THREE.Plane(new THREE.Vector3(0, 0, 1), 0), hit)) {
    sendGoal(hit.x, hit.y, null);
  }
});

// ============================================================================ goal

function toast(text, bad = false) {
  const t = $("toast");
  t.textContent = text;
  t.style.color = bad ? css("--critical") : "";
  t.hidden = false;
  clearTimeout(toast.timer);
  toast.timer = setTimeout(() => { t.hidden = true; }, 3000);
}
async function post(url, body) {
  const r = await fetch(url, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
  return r.json();
}
async function sendGoal(x, y, z) {
  const zz = z === null || z === "" || Number.isNaN(z) ? info.defaultGoalZ : Number(z);
  try {
    const res = await post("/api/goal", { x, y, z: zz });
    toast(res.ok ? `Goal sent: (${fmt(x)}, ${fmt(y)}, ${fmt(zz)})` : `Goal not sent: ${res.reason}`, !res.ok);
  } catch (e) {
    toast(`Goal not sent: ${e}`, true);
  }
}
$("goal-form").addEventListener("submit", (e) => {
  e.preventDefault();
  const f = new FormData(e.target);
  sendGoal(Number(f.get("x")), Number(f.get("y")), f.get("z") === "" ? null : Number(f.get("z")));
});

// ============================================================================ tiles

function tiles(el, items) {
  el.innerHTML = items.map(([label, value, unit, sub]) =>
    `<div class="tile"><div class="label">${label}</div><div class="value">${value}<span class="unit">${unit || ""}</span></div>${sub ? `<div class="sub">${sub}</div>` : ""}</div>`).join("");
}

const statusLabels = {
  YAWING: "Yawing", TRAVELING: "Traveling", GOAL_SEEN: "Goal seen",
  GOAL_REACHED: "Goal reached", HOVER_AVOIDING: "Hover avoiding", UNKNOWN: "Waiting for SANDO",
};
function updateState(s) {
  const st = s.connected ? s.status : "OFFLINE";
  $("status").className = `status status-${st.toLowerCase()}`;
  $("status-text").textContent = s.connected ? (statusLabels[st] || st) : "No status from SANDO";
  if (s.pos) {
    vehicle.position.set(...s.pos);
    vehicle.quaternion.set(...s.quat);
    if (!viewInitialized) { viewInitialized = true; lastVehiclePos.copy(vehicle.position); resetView(); }
  }
  pointMeshes.goal.visible = !!s.goal;
  if (s.goal) pointMeshes.goal.position.set(...s.goal);
  for (const k of ["A", "G", "E"]) {
    const p = s.points?.[k];
    pointMeshes[k].visible = !!p;
    if (p) pointMeshes[k].position.set(...p);
  }
  tiles($("vehicle-tiles"), [
    ["Speed", fmt(s.speed, 2), " m/s"],
    ["Distance to goal", fmt(s.goalDist, 1), " m"],
    ["Altitude", s.pos ? fmt(s.pos[2], 2) : "—", " m"],
    ["Position x, y", s.pos ? `${fmt(s.pos[0])}, ${fmt(s.pos[1])}` : "—", ""],
    ["Goal", s.goal ? `${fmt(s.goal[0])}, ${fmt(s.goal[1])}` : "—", "", s.goal ? `z ${fmt(s.goal[2])}` : ""],
    ["Vertical speed", s.vel ? fmt(s.vel[2], 2) : "—", " m/s"],
  ]);
}

// ============================================================================ chart

const SERIES = [
  { key: "total", label: "Total", color: "--series-1" },
  { key: "local", label: "Local trajectory", color: "--series-2" },
  { key: "global", label: "Global path", color: "--series-3" },
  { key: "decomp", label: "Corridors", color: "--series-4" },
];
$("legend").innerHTML = SERIES.map((s) => `<span><i class="swatch" style="background:var(${s.color})"></i>${s.label}</span>`).join("")
  + `<span><i class="swatch" style="background:var(--critical);width:2px;height:10px"></i>Failed replan</span>`;

function niceMax(v) {
  if (!(v > 0)) return 10;
  const p = 10 ** Math.floor(Math.log10(v));
  for (const m of [1, 2, 2.5, 5, 10]) if (m * p >= v) return m * p;
  return 10 * p;
}

const chart = {
  canvas: $("chart"),
  data: null,
  hoverX: null,
  pad: { l: 40, r: 8, t: 8, b: 30 },
  draw() {
    const c = this.canvas, d = this.data;
    const dpr = window.devicePixelRatio || 1;
    const w = c.clientWidth, h = c.clientHeight;
    c.width = w * dpr; c.height = h * dpr;
    const g = c.getContext("2d");
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, w, h);
    const { l, r, t, b } = this.pad;
    const pw = w - l - r, ph = h - t - b;
    const xOf = (s) => l + ((s + 60) / 60) * pw;
    let ymax = 10;
    if (d && d.t.length) {
      const vals = d.total.filter((v) => v > 0).sort((a, b) => a - b);
      if (vals.length) ymax = niceMax(vals[Math.floor(0.98 * (vals.length - 1))] * 1.1);
    }
    this.ymax = ymax;
    const yOf = (v) => t + ph - (Math.min(v, ymax) / ymax) * ph;

    g.font = "11px system-ui, sans-serif";
    g.fillStyle = css("--muted");
    g.strokeStyle = css("--grid");
    g.lineWidth = 1;
    g.textAlign = "right"; g.textBaseline = "middle";
    for (let i = 0; i <= 4; i++) {
      const v = (ymax * i) / 4, y = Math.round(yOf(v)) + 0.5;
      g.beginPath(); g.moveTo(l, y); g.lineTo(l + pw, y); g.stroke();
      g.fillText(v % 1 ? v.toFixed(1) : String(v), l - 6, y);
    }
    g.textAlign = "center"; g.textBaseline = "top";
    for (const s of [-60, -45, -30, -15, 0]) g.fillText(s === 0 ? "now" : `${s} s`, xOf(s), t + ph + 14);
    if (!d || !d.t.length) {
      g.fillText("No replans yet", l + pw / 2, t + ph / 2);
      return;
    }
    // Failed replans: ticks just above the baseline
    g.strokeStyle = css("--critical");
    g.lineWidth = 1.5;
    for (let i = 0; i < d.t.length; i++) {
      if (d.ok[i]) continue;
      const x = xOf(d.t[i]);
      g.beginPath(); g.moveTo(x, t + ph + 2); g.lineTo(x, t + ph + 10); g.stroke();
    }
    g.lineWidth = 2; g.lineJoin = "round";
    for (const s of [...SERIES].reverse()) {
      g.strokeStyle = css(s.color);
      g.beginPath();
      let pen = false;
      for (let i = 0; i < d.t.length; i++) {
        const v = d[s.key][i];
        if (!d.ok[i] && s.key !== "total") { pen = false; continue; }
        const x = xOf(d.t[i]), y = yOf(v);
        if (pen) g.lineTo(x, y); else { g.moveTo(x, y); pen = true; }
      }
      g.stroke();
    }
    g.strokeStyle = css("--axis"); g.lineWidth = 1;
    g.beginPath(); g.moveTo(l, t + ph + 0.5); g.lineTo(l + pw, t + ph + 0.5); g.stroke();

    if (this.hoverX !== null) {
      const sec = ((this.hoverX - l) / pw) * 60 - 60;
      let best = 0;
      for (let i = 1; i < d.t.length; i++) if (Math.abs(d.t[i] - sec) < Math.abs(d.t[best] - sec)) best = i;
      const x = xOf(d.t[best]);
      g.strokeStyle = css("--axis");
      g.beginPath(); g.moveTo(x + 0.5, t); g.lineTo(x + 0.5, t + ph); g.stroke();
      for (const s of SERIES) {
        if (!d.ok[best] && s.key !== "total") continue;
        g.fillStyle = css(s.color);
        g.strokeStyle = css("--surface");
        g.lineWidth = 2;
        g.beginPath(); g.arc(x, yOf(d[s.key][best]), 4, 0, 2 * Math.PI); g.fill(); g.stroke();
      }
      const tip = $("tooltip");
      tip.innerHTML = `<div><b>${fmt(-d.t[best], 1)} s ago</b> · ${d.ok[best] ? "succeeded" : "<span style='color:var(--critical)'>failed</span>"}</div>`
        + SERIES.map((s) => `<div><i class="swatch" style="background:var(${s.color})"></i> <span class="k">${s.label}</span> ${fmt(d[s.key][best], 1)} ms</div>`).join("");
      tip.hidden = false;
      const tw = tip.offsetWidth;
      tip.style.left = `${x + 12 + tw > w ? x - tw - 12 : x + 12}px`;
      tip.style.top = `${t}px`;
    }
  },
};
chart.canvas.addEventListener("mousemove", (e) => {
  chart.hoverX = e.offsetX;
  chart.draw();
});
chart.canvas.addEventListener("mouseleave", () => {
  chart.hoverX = null;
  $("tooltip").hidden = true;
  chart.draw();
});
new ResizeObserver(() => chart.draw()).observe(chart.canvas);

const timingLabels = {
  total_replanning_ms: "Total", global_planning_ms: "Global path", cvx_decomp_ms: "Corridors",
  local_traj_ms: "Local trajectory", safe_paths_ms: "Safe paths", safety_check_ms: "Safety check",
  yaw_sequence_ms: "Yaw sequence",
};
let perfCache = null;
function updatePerf(p) {
  perfCache = p;
  const tot = p.stats.total_replanning_ms, loc = p.stats.local_traj_ms;
  tiles($("perf-tiles"), [
    ["Replans", fmt(p.replanHz, 0), " Hz", `${fmt(p.successHz, 0)} Hz succeeded`],
    ["Success", fmt(p.successPct, 0), " %", `${p.samples} replans`],
    ["Total time", fmt(tot?.median, 1), " ms", tot ? `p95 ${fmt(tot.p95, 1)} · max ${fmt(tot.max, 0)}` : ""],
    ["Local trajectory", fmt(loc?.median, 1), " ms", loc ? `p95 ${fmt(loc.p95, 1)}` : ""],
    ["Time factor", fmt(p.lastFactor, 2), "", "last success"],
    ["Goals reached", p.goalsReached, ""],
  ]);
  $("timing-table").innerHTML = `<tr><th>Step (successful replans)</th><th class="num">Median</th><th class="num">p95</th><th class="num">Max</th></tr>`
    + Object.entries(p.stats).map(([k, s]) => `<tr><td>${timingLabels[k] || k}</td><td class="num">${fmt(s?.median, 2)}</td><td class="num">${fmt(s?.p95, 2)}</td><td class="num">${fmt(s?.max, 1)}</td></tr>`).join("");
  chart.data = p.series;
  chart.draw();
}

// ============================================================================ parameters

let params = [];
function paramMsg(text, bad = false) {
  const m = $("param-msg");
  m.textContent = text;
  m.className = `small ${bad ? "bad" : "muted"}`;
}
async function loadParams() {
  try {
    const r = await fetch("/api/params");
    const d = await r.json();
    if (!r.ok) throw new Error(d.error);
    params = d.params;
    $("param-node").textContent = d.node;
    paramMsg(d.readOnly ? "Read-only: the server was started with read_only:=true" : "");
    renderParams();
  } catch (e) {
    paramMsg(`Parameters unavailable: ${e.message}`, true);
  }
}
function renderParams() {
  const q = $("param-filter").value.trim().toLowerCase();
  const show = params.filter((p) => !q || p.name.toLowerCase().includes(q));
  const editable = !info.readOnly;
  const row = (p) => {
    const tr = document.createElement("tr");
    const name = document.createElement("td");
    name.textContent = p.name;
    const val = document.createElement("td");
    if (p.live && editable) {
      let input;
      if (p.type === "bool") {
        input = document.createElement("input");
        input.type = "checkbox";
        input.checked = p.value;
        input.addEventListener("change", () => setParam(p, input.checked, tr));
      } else {
        input = document.createElement("input");
        input.type = p.type === "str" ? "text" : "number";
        if (p.type === "float") input.step = "any";
        input.value = p.value;
        input.addEventListener("keydown", (e) => { if (e.key === "Enter") setParam(p, input.value, tr); });
        input.addEventListener("change", () => setParam(p, input.value, tr));
      }
      input.setAttribute("aria-label", p.name);
      val.appendChild(input);
    } else {
      val.className = "ro";
      val.textContent = Array.isArray(p.value) ? `[${p.value.join(", ")}]` : String(p.value);
    }
    tr.append(name, val);
    return tr;
  };
  $("params-live").replaceChildren(...show.filter((p) => p.live).map(row));
  $("params-static").replaceChildren(...show.filter((p) => !p.live).map(row));
}
async function setParam(p, value, tr) {
  try {
    const res = await post("/api/params", { name: p.name, value });
    if (res.ok) {
      p.value = p.type === "bool" ? value : Number.isNaN(Number(value)) ? value : Number(value);
      tr.classList.remove("err");
      paramMsg(`${p.name} set to ${value}`);
    } else {
      tr.classList.add("err");
      paramMsg(`${p.name} not changed: ${res.reason}`, true);
    }
  } catch (e) {
    paramMsg(`${p.name} not changed: ${e}`, true);
  }
}
$("param-filter").addEventListener("input", renderParams);
$("param-refresh").addEventListener("click", loadParams);

// ============================================================================ events

const eventsEl = $("events");
function addEvents(list) {
  for (const e of list) {
    const li = document.createElement("li");
    const t = document.createElement("time");
    t.textContent = new Date(e.t * 1000).toLocaleTimeString();
    li.append(t, document.createTextNode(e.text));
    eventsEl.prepend(li);
    if (e.kind === "param") loadParams();
  }
  while (eventsEl.children.length > 200) eventsEl.lastChild.remove();
}

// ============================================================================ stream

function connect() {
  const es = new EventSource("/events");
  const conn = $("conn");
  let lastMsg = Date.now();
  es.onopen = () => { conn.textContent = "live"; conn.className = "conn ok"; };
  es.onerror = () => { conn.textContent = "disconnected, retrying…"; conn.className = "conn bad"; };
  const on = (kind, fn) => es.addEventListener(kind, (e) => {
    lastMsg = Date.now();
    try { fn(JSON.parse(e.data)); } catch (err) { console.error(kind, err); }
  });
  on("state", updateState);
  on("scene", (s) => { updateMarkers(s.markers); updateTrail(s.trail); updateObstacles(s.obstacles); });
  on("polys", updatePolys);
  on("clouds", updateClouds);
  on("perf", updatePerf);
  on("events", addEvents);
  setInterval(() => {
    if (Date.now() - lastMsg > 3000 && es.readyState === EventSource.OPEN) {
      conn.textContent = "connected, no data"; conn.className = "conn bad";
    } else if (es.readyState === EventSource.OPEN) {
      conn.textContent = "live"; conn.className = "conn ok";
    }
  }, 1000);
}

(async function start() {
  try {
    info = await (await fetch("/api/info")).json();
    $("agent").textContent = `/${info.namespace}`;
    document.title = `SANDO Monitor · ${info.namespace}`;
  } catch (e) { /* server will report via the stream */ }
  connect();
  loadParams();
})();
