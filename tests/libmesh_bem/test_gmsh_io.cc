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
 * @file test_gmsh_io.cc
 * @brief Verify the Gmsh MSH 4.1 binary reader against libMesh's ASCII reader.
 */

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <libmesh/elem.h>
#include <libmesh/replicated_mesh.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "libmesh_bem/gmsh_io.h"
#include "test_common.h"

using namespace HierBEM::LibMeshBEM;
using libMesh::Point;
using libMesh::Real;

namespace
{
  /**
   * Collect the vertex coordinates of all triangles in a canonical order, so
   * that meshes with different node and element numberings can be compared.
   */
  std::vector<std::vector<Real>>
  canonical_triangles(const libMesh::MeshBase &mesh)
  {
    std::vector<std::vector<Real>> tris;
    for (const libMesh::Elem *elem : mesh.active_element_ptr_range())
      {
        if (elem->dim() != 2)
          continue;
        std::vector<std::vector<Real>> v;
        for (unsigned int i = 0; i < elem->n_vertices(); ++i)
          v.push_back(
            {elem->point(i)(0), elem->point(i)(1), elem->point(i)(2)});
        std::sort(v.begin(), v.end());
        std::vector<Real> flat;
        for (const auto &p : v)
          flat.insert(flat.end(), p.begin(), p.end());
        tris.push_back(flat);
      }
    std::sort(tris.begin(), tris.end());
    return tris;
  }
} // namespace


TEST_CASE("Read Gmsh MSH 4.1 binary file", "[libmesh_bem]")
{
  libMesh::ReplicatedMesh binary_mesh(
    HierBEM::LibMeshBEM::testing::libmesh_init->comm(), 3);
  read_gmsh(binary_mesh, TEST_DIR "/unit_square_coarse_binary.msh");

  libMesh::ReplicatedMesh ascii_mesh(
    HierBEM::LibMeshBEM::testing::libmesh_init->comm(), 3);
  read_gmsh(ascii_mesh, TEST_DIR "/unit_square_coarse_ascii.msh");

  REQUIRE(binary_mesh.mesh_dimension() == 2);
  REQUIRE(binary_mesh.n_elem() > 0);
  REQUIRE(binary_mesh.n_elem() == ascii_mesh.n_elem());
  REQUIRE(binary_mesh.n_nodes() == ascii_mesh.n_nodes());

  Real area = 0.;
  for (const libMesh::Elem *elem : binary_mesh.active_element_ptr_range())
    {
      REQUIRE(elem->type() == libMesh::TRI3);
      REQUIRE(elem->subdomain_id() == 1);
      area += elem->volume();
    }
  REQUIRE(area == Catch::Approx(1.0).epsilon(1e-12));
  REQUIRE(binary_mesh.subdomain_name(1) == "plate");

  // The ASCII file stores the coordinates with 16 significant digits only.
  const auto binary_tris = canonical_triangles(binary_mesh);
  const auto ascii_tris  = canonical_triangles(ascii_mesh);
  REQUIRE(binary_tris.size() == ascii_tris.size());
  for (std::size_t i = 0; i < binary_tris.size(); ++i)
    for (std::size_t j = 0; j < binary_tris[i].size(); ++j)
      REQUIRE(std::abs(binary_tris[i][j] - ascii_tris[i][j]) < 1e-14);
}


TEST_CASE("Read the example mesh", "[libmesh_bem]")
{
  libMesh::ReplicatedMesh mesh(
    HierBEM::LibMeshBEM::testing::libmesh_init->comm(), 3);
  read_gmsh41_binary(mesh, EXAMPLE_MESH);

  Real area = 0.;
  for (const libMesh::Elem *elem : mesh.active_element_ptr_range())
    {
      REQUIRE(elem->type() == libMesh::TRI3);
      for (unsigned int i = 0; i < 3; ++i)
        REQUIRE(elem->point(i)(2) == 0.);
      area += elem->volume();
    }
  REQUIRE(area == Catch::Approx(1.0).epsilon(1e-12));
}
