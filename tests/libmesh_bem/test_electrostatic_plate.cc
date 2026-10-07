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
 * @file test_electrostatic_plate.cc
 * @brief Capacitance of the unit square plate computed with P0 Galerkin BEM.
 */

#include <libmesh/replicated_mesh.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "libmesh_bem/galerkin_assembly.h"
#include "libmesh_bem/gmsh_io.h"
#include "libmesh_bem/laplace_kernels.h"
#include "test_common.h"

using namespace HierBEM::LibMeshBEM;
using libMesh::Real;

TEST_CASE("Capacitance of the unit square plate", "[libmesh_bem]")
{
  libMesh::ReplicatedMesh mesh(
    HierBEM::LibMeshBEM::testing::libmesh_init->comm(), 3);
  read_gmsh(mesh, EXAMPLE_MESH);

  const auto elems = collect_p0_elements(mesh);
  const auto areas = element_areas(elems);

  const SauterTriangleQuadrature quad;
  libMesh::DenseMatrix<Real>     V =
    assemble_p0_matrix(elems, LaplaceSingleLayerKernel(), quad, true);
  libMesh::DenseVector<Real> b =
    assemble_p0_rhs(elems, [](const libMesh::Point &) { return 1.; }, 2);

  libMesh::DenseVector<Real> s;
  V.lu_solve(b, s);

  // Normalized capacitance \f$C / (4\pi\varepsilon_0)\f$ for unit potential.
  Real q = 0.;
  for (unsigned int i = 0; i < s.size(); ++i)
    {
      // The charge density must be positive everywhere on the plate.
      REQUIRE(s(i) > 0.);
      q += s(i) * areas[i];
    }
  const Real normalized_capacitance = q / (4. * libMesh::pi);

  // Reference value 0.3667874; the coarse P0 discretization underestimates
  // it slightly.
  REQUIRE(normalized_capacitance == Catch::Approx(0.3667874).epsilon(0.02));
  REQUIRE(normalized_capacitance < 0.3667874);
}
