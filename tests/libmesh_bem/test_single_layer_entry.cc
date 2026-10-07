// Copyright (C) 2026 HierBEM contributors
//
// This file is part of the HierBEM library.
//
// HierBEM is free software: you can use it, redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version. The full text of the license can be found in the
// file LICENSE at the top level directory of HierBEM.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <libmesh/equation_systems.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/system.h>

#include "libmesh_bem/galerkin_assembly.h"
#include "libmesh_bem/gmsh_io.h"
#include "libmesh_bem/laplace_kernels.h"
#include "libmesh_bem/single_layer_entry.h"
#include "test_common.h"

using namespace HierBEM::LibMeshBEM;
using libMesh::Real;

namespace
{
  void
  check_entries(const SingleLayerEntryEvaluator         &evaluator,
                const libMesh::DenseMatrix<libMesh::Real> &V)
  {
    REQUIRE(evaluator.n_dofs() == V.m());
    for (unsigned int i = 0; i < V.m(); ++i)
      for (unsigned int j = 0; j < V.n(); ++j)
        {
          INFO("Entry " << i << ", " << j);
          REQUIRE(evaluator.entry(i, j) ==
                  Catch::Approx(V(i, j)).epsilon(1e-12));

          // ButterflyPACK passes one-based indices and a void* user pointer.
          int    m = i + 1, n = j + 1;
          double val = 0.;
          butterflypack_single_layer_element(
            &m, &n, &val, const_cast<SingleLayerEntryEvaluator *>(&evaluator));
          REQUIRE(val == evaluator.entry(i, j));
        }
  }
} // namespace

TEST_CASE("Single layer entries match dense P0 assembly", "[libmesh_bem]")
{
  libMesh::ReplicatedMesh mesh(testing::libmesh_init->comm(), 3);
  read_gmsh(mesh, TEST_DIR "/unit_square_coarse_ascii.msh");
  const auto                     elems = collect_p0_elements(mesh);
  const SauterTriangleQuadrature quad;
  const auto                     V =
    assemble_p0_matrix(elems, LaplaceSingleLayerKernel(), quad, false);
  check_entries(SingleLayerEntryEvaluator(elems, quad), V);
}

TEST_CASE("Single layer entries match dense P1 and P2 assembly",
          "[libmesh_bem]")
{
  for (const auto order : {libMesh::FIRST, libMesh::SECOND})
    {
      INFO("Order " << order);
      libMesh::ReplicatedMesh mesh(testing::libmesh_init->comm(), 3);
      read_gmsh(mesh, TEST_DIR "/unit_square_coarse_ascii.msh");
      if (order == libMesh::SECOND)
        mesh.all_second_order();
      libMesh::EquationSystems systems(mesh);
      auto &system = systems.add_system<libMesh::System>("charge");
      system.add_variable("s", order, libMesh::LAGRANGE);
      systems.init();
      const LagrangeTriangleSpace    space(mesh, system.get_dof_map());
      const SauterTriangleQuadrature quad;
      const auto                     V =
        assemble_lagrange_matrix(space, LaplaceSingleLayerKernel(), quad, false);
      check_entries(SingleLayerEntryEvaluator(space, quad), V);
    }
}
