# Copyright (C) 2026 HierBEM contributors
#
# This file is part of the HierBEM library.
#
# HierBEM is free software: you can use it, redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version. The full text of the license can be found in the
# file LICENSE at the top level directory of HierBEM.
"""Generate a triangular surface mesh of the unit square [0,1]x[0,1]x{0}.

The mesh is written in the Gmsh MSH 4.1 *binary* format.

Usage:
    python3 generate_unit_square_mesh.py [output.msh] [mesh_size]
"""

import sys

import gmsh


def main():
    output = sys.argv[1] if len(sys.argv) > 1 else "unit_square.msh"
    mesh_size = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("unit_square")

    p = [
        gmsh.model.geo.addPoint(0, 0, 0, mesh_size),
        gmsh.model.geo.addPoint(1, 0, 0, mesh_size),
        gmsh.model.geo.addPoint(1, 1, 0, mesh_size),
        gmsh.model.geo.addPoint(0, 1, 0, mesh_size),
    ]
    lines = [gmsh.model.geo.addLine(p[i], p[(i + 1) % 4]) for i in range(4)]
    loop = gmsh.model.geo.addCurveLoop(lines)
    surface = gmsh.model.geo.addPlaneSurface([loop])
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(2, [surface], 1, name="plate")

    # Linear triangles only.
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.RecombineAll", 0)
    gmsh.model.mesh.generate(2)

    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(output)
    gmsh.finalize()


if __name__ == "__main__":
    main()
