#!/usr/bin/env python3
# /* ----------------------------------------------------------------------------
#  * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

"""Convert a Gazebo Classic .world file into a Gazebo Sim (gz) world for worlds/gz/.

- Drops the Classic-only plugins (gazebo_ros_state, disable_physics), <gui>, <audio> and the
  Classic <physics> block, and adds the gz systems SANDO needs (physics, user commands for
  set_pose/create, scene broadcaster, sensors).
- Applies the model poses saved in <state> (Classic uses them over the <model> poses).
- Replaces Gazebo/<Color> material scripts (OGRE 1 only) with plain colors.
- Replaces <include> of model://sun and model://ground_plane with equivalent inline elements.

Usage:
    python3 scripts/convert_world_to_gz.py worlds/easy_forest.world worlds/gz/easy_forest.sdf
"""

import argparse
import sys
import xml.etree.ElementTree as ET

# Colors of the Classic Gazebo/<Name> materials (media/materials/scripts/gazebo.material)
GAZEBO_COLORS = {
    "Gazebo/Grey": ((0.3, 0.3, 0.3, 1), (0.7, 0.7, 0.7, 1)),
    "Gazebo/DarkGrey": ((0.175, 0.175, 0.175, 1), (0.175, 0.175, 0.175, 1)),
    "Gazebo/White": ((1, 1, 1, 1), (1, 1, 1, 1)),
    "Gazebo/Black": ((0, 0, 0, 1), (0, 0, 0, 1)),
    "Gazebo/Red": ((1, 0, 0, 1), (1, 0, 0, 1)),
    "Gazebo/Green": ((0, 1, 0, 1), (0, 1, 0, 1)),
    "Gazebo/Blue": ((0, 0, 1, 1), (0, 0, 1, 1)),
    "Gazebo/Yellow": ((1, 1, 0, 1), (1, 1, 0, 1)),
    "Gazebo/Orange": ((1, 0.5088, 0.0468, 1), (1, 0.5088, 0.0468, 1)),
    "Gazebo/Wood": ((0.5, 0.35, 0.2, 1), (0.5, 0.35, 0.2, 1)),
}

GZ_SYSTEMS = [
    # Applies the pose commands sent by fake_sim (set_pose)
    ("gz-sim-physics-system", "gz::sim::systems::Physics", None),
    # set_pose and create services
    ("gz-sim-user-commands-system", "gz::sim::systems::UserCommands", None),
    ("gz-sim-scene-broadcaster-system", "gz::sim::systems::SceneBroadcaster", None),
    # Rendering-based sensors (D435 depth/color cameras)
    ("gz-sim-sensors-system", "gz::sim::systems::Sensors", ("render_engine", "ogre2")),
]

SUN_LIGHT = """<light name='sun' type='directional'>
  <cast_shadows>1</cast_shadows>
  <pose>0 0 10 0 0 0</pose>
  <diffuse>0.8 0.8 0.8 1</diffuse>
  <specular>0.2 0.2 0.2 1</specular>
  <attenuation>
    <range>1000</range>
    <constant>0.9</constant>
    <linear>0.01</linear>
    <quadratic>0.001</quadratic>
  </attenuation>
  <direction>-0.5 0.1 -0.9</direction>
</light>"""

GROUND_PLANE = """<model name='ground_plane'>
  <static>true</static>
  <link name='link'>
    <collision name='collision'>
      <geometry><plane><normal>0 0 1</normal><size>100 100</size></plane></geometry>
    </collision>
    <visual name='visual'>
      <cast_shadows>false</cast_shadows>
      <geometry><plane><normal>0 0 1</normal><size>100 100</size></plane></geometry>
      <material>
        <ambient>0.8 0.8 0.8 1</ambient>
        <diffuse>0.8 0.8 0.8 1</diffuse>
        <specular>0.8 0.8 0.8 1</specular>
      </material>
    </visual>
  </link>
</model>"""


def _fmt(color):
    return " ".join(f"{c:g}" for c in color)


def convert(root, source_name):
    world = root.find("world")
    if world is None:
        raise ValueError("no <world> element")

    # Model poses saved in <state> win over the <model> poses in Gazebo Classic
    state_poses = {}
    for state in world.findall("state"):
        for model in state.findall("model"):
            pose = model.find("pose")
            if pose is not None and pose.text:
                state_poses[model.get("name")] = pose.text.strip()

    for tag in ("plugin", "gui", "state", "audio", "physics"):
        for el in world.findall(tag):
            world.remove(el)

    # model://sun and model://ground_plane have no local copy in gz: inline what they contain
    for inc in list(world):
        if inc.tag != "include":
            continue
        uri = (inc.findtext("uri") or "").strip()
        i = list(world).index(inc)
        if uri == "model://sun":
            world.remove(inc)
            if world.find("light[@name='sun']") is None:
                world.insert(i, ET.fromstring(SUN_LIGHT))
        elif uri == "model://ground_plane":
            world.remove(inc)
            ground = ET.fromstring(GROUND_PLANE)
            for tag in ("name", "pose"):
                if inc.find(tag) is not None:
                    if tag == "name":
                        ground.set("name", inc.findtext("name").strip())
                    else:
                        ground.insert(0, inc.find("pose"))
            world.insert(i, ground)
        else:
            print(f"[WARN] keeping <include> {uri}: it must be on GZ_SIM_RESOURCE_PATH",
                  file=sys.stderr)

    # gz systems and a 10 ms step (Classic ran 1 ms steps with physics disabled)
    head = []
    physics = ET.Element("physics", name="10ms", type="ignored")
    ET.SubElement(physics, "max_step_size").text = "0.01"
    ET.SubElement(physics, "real_time_factor").text = "1"
    head.append(physics)
    for filename, name, param in GZ_SYSTEMS:
        plugin = ET.Element("plugin", filename=filename, name=name)
        if param:
            ET.SubElement(plugin, param[0]).text = param[1]
        head.append(plugin)
    for i, el in enumerate(head):
        world.insert(i, el)

    for model in world.iter("model"):
        name = model.get("name")
        if name in state_poses:
            pose = model.find("pose")
            if pose is None:
                pose = ET.SubElement(model, "pose")
            pose.text = state_poses[name]

    unknown = set()
    for material in world.iter("material"):
        script = material.find("script")
        if script is None:
            continue
        mat_name = (script.findtext("name") or "").strip()
        material.remove(script)
        if mat_name not in GAZEBO_COLORS:
            unknown.add(mat_name)
            continue
        ambient, diffuse = GAZEBO_COLORS[mat_name]
        for tag, value in (("ambient", ambient), ("diffuse", diffuse)):
            el = material.find(tag)
            if el is None:
                el = ET.SubElement(material, tag)
            el.text = _fmt(value)
        if material.find("specular") is None:
            ET.SubElement(material, "specular").text = "0.1 0.1 0.1 1"
    for mat_name in sorted(unknown):
        print(f"[WARN] unknown material script '{mat_name}': left without a color", file=sys.stderr)

    root.set("version", "1.9")
    root.insert(0, ET.Comment(f" Gazebo Sim (gz) version of worlds/{source_name}, "
                              "generated by scripts/convert_world_to_gz.py "))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="Gazebo Classic .world file")
    parser.add_argument("output", help="Gazebo Sim .sdf file to write")
    args = parser.parse_args()

    tree = ET.parse(args.input)
    convert(tree.getroot(), args.input.rsplit("/", 1)[-1])
    ET.indent(tree, space="  ")
    tree.write(args.output, encoding="unicode", xml_declaration=True)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
