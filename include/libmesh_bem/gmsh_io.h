// Copyright (C) 2026 HierBEM contributors
//
// This file is part of the HierBEM library.
//
// HierBEM is free software: you can use it, redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version. The full text of the license can be found in the
// file LICENSE at the top level directory of HierBEM.

/**
 * @file gmsh_io.h
 * @brief Reader for Gmsh MSH 4.1 files (binary and ASCII) into libMesh
 * meshes.
 *
 * libMesh's own @p GmshIO only supports the ASCII variant of the MSH format,
 * while Gmsh writes binary files when @p Mesh.Binary=1. This reader supports
 * the binary MSH 4.1 format and delegates ASCII files to libMesh's @p GmshIO.
 */

#ifndef HIERBEM_INCLUDE_LIBMESH_BEM_GMSH_IO_H_
#define HIERBEM_INCLUDE_LIBMESH_BEM_GMSH_IO_H_

#include <libmesh/mesh_base.h>

#include <string>

namespace HierBEM
{
  namespace LibMeshBEM
  {
    /**
     * Read the surface (2D) elements of a binary Gmsh MSH 4.1 file into
     * @p mesh. Points and curves are ignored. Each element gets the physical
     * tag of its geometric entity as subdomain id (or the entity tag if no
     * physical group is defined), and the physical names are set as subdomain
     * names. Supported element types are 3-node triangles, 6-node triangles
     * and 4-node quadrilaterals.
     *
     * The mesh is cleared before reading and @p prepare_for_use() is called
     * at the end.
     */
    void
    read_gmsh41_binary(libMesh::MeshBase &mesh, const std::string &filename);

    /**
     * Read a Gmsh file. Binary MSH 4.1 files are read with
     * @p read_gmsh41_binary, ASCII files with libMesh's @p GmshIO.
     */
    void
    read_gmsh(libMesh::MeshBase &mesh, const std::string &filename);
  } // namespace LibMeshBEM
} // namespace HierBEM

#endif // HIERBEM_INCLUDE_LIBMESH_BEM_GMSH_IO_H_
