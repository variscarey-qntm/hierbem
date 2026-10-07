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

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

#include "libmesh_bem/galerkin_assembly.h"
#include "libmesh_bem/gmsh_io.h"
#include "libmesh_bem/laplace_kernels.h"
#include "test_common.h"

using namespace HierBEM::LibMeshBEM;
using libMesh::Point;
using libMesh::Real;

TEST_CASE("Lagrange triangle DOFs and exact basis integrals", "[libmesh_bem]")
{
  for (const auto &[order, quadratic_mesh] :
       {std::make_pair(libMesh::FIRST, false),
        std::make_pair(libMesh::FIRST, true),
        std::make_pair(libMesh::SECOND, true)})
    {
      INFO("Order " << order);
      libMesh::ReplicatedMesh mesh(testing::libmesh_init->comm(), 3);
      for (const Point p :
           {Point(0., 0.), Point(1., 0.), Point(1., 1.), Point(0., 1.)})
        mesh.add_point(p);
      for (const auto nodes : {std::array<unsigned int, 3>{{0, 1, 2}},
                               std::array<unsigned int, 3>{{3, 2, 0}}})
        {
          auto elem = libMesh::Elem::build(libMesh::TRI3);
          for (unsigned int i = 0; i < 3; ++i)
            elem->set_node(i) = mesh.node_ptr(nodes[i]);
          mesh.add_elem(std::move(elem));
        }
      mesh.prepare_for_use();
      if (quadratic_mesh)
        mesh.all_second_order();
      libMesh::EquationSystems systems(mesh);
      auto &system = systems.add_system<libMesh::System>("charge");
      system.add_variable("s", order, libMesh::LAGRANGE);
      systems.init();
      const LagrangeTriangleSpace space(mesh, system.get_dof_map());
      REQUIRE(space.n_dofs() == (order == libMesh::FIRST ? 4 : 9));
      const auto  &d0     = space.dof_indices(0);
      const auto  &d1     = space.dof_indices(1);
      unsigned int shared = 0;
      for (const auto d : d0)
        shared += std::count(d1.begin(), d1.end(), d);
      REQUIRE(shared == (order == libMesh::FIRST ? 2 : 3));

      const auto b =
        assemble_lagrange_rhs(space, [](const Point &) { return 1.; }, 3);
      std::vector<Real> exact(space.n_dofs(), 0.);
      for (unsigned int e = 0; e < 2; ++e)
        for (unsigned int a = 0; a < space.dof_indices(e).size(); ++a)
          exact[space.dof_indices(e)[a]] +=
            (order == libMesh::FIRST || a >= 3) ? 1. / 6. : 0.;
      for (unsigned int i = 0; i < b.size(); ++i)
        REQUIRE(b(i) == Catch::Approx(exact[i]).margin(1e-13));
    }
}

TEST_CASE("Lagrange assembly respects permutations and partition of unity",
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
      const auto f  = [](const Point &x) { return 1. + x(0); };
      const auto g  = [](const Point &y) { return 2. - y(1); };
      const auto bf = assemble_lagrange_rhs(space, f, 4);
      const auto bg = assemble_lagrange_rhs(space, g, 4);
      const auto A  = assemble_lagrange_matrix(
        space,
        [&](const Point &x, const Point &y, const Point &, const Point &) {
          return f(x) * g(y);
        },
        quad,
        false);
      for (unsigned int i = 0; i < A.m(); ++i)
        for (unsigned int j = 0; j < A.n(); ++j)
          REQUIRE(A(i, j) == Catch::Approx(bf(i) * bg(j)).margin(1e-12));

      const auto V =
        assemble_lagrange_matrix(space, LaplaceSingleLayerKernel(), quad, true);
      const auto full = assemble_lagrange_matrix(space,
                                                 LaplaceSingleLayerKernel(),
                                                 quad,
                                                 false);
      const auto p0   = assemble_p0_matrix(space.elements(),
                                         LaplaceSingleLayerKernel(),
                                         quad,
                                         true);
      Real       sum = 0., sum_p0 = 0.;
      for (unsigned int i = 0; i < V.m(); ++i)
        for (unsigned int j = 0; j < V.n(); ++j)
          {
            REQUIRE(V(i, j) == Catch::Approx(full(i, j)).margin(1e-12));
            REQUIRE(V(i, j) == Catch::Approx(V(j, i)).margin(1e-12));
            sum += V(i, j);
          }
      for (unsigned int i = 0; i < p0.m(); ++i)
        for (unsigned int j = 0; j < p0.n(); ++j)
          sum_p0 += p0(i, j);
      REQUIRE(sum == Catch::Approx(sum_p0).epsilon(1e-12));
    }
}

TEST_CASE("Dense P1 and P2 square plate capacitance", "[libmesh_bem]")
{
  for (const auto order : {libMesh::FIRST, libMesh::SECOND})
    {
      INFO("Order " << order);
      libMesh::ReplicatedMesh mesh(testing::libmesh_init->comm(), 3);
      read_gmsh(mesh, EXAMPLE_MESH);
      if (order == libMesh::SECOND)
        mesh.all_second_order();
      libMesh::EquationSystems systems(mesh);
      auto &system = systems.add_system<libMesh::System>("charge");
      system.add_variable("s", order, libMesh::LAGRANGE);
      systems.init();
      const LagrangeTriangleSpace    space(mesh, system.get_dof_map());
      const SauterTriangleQuadrature quad;
      auto                           V =
        assemble_lagrange_matrix(space, LaplaceSingleLayerKernel(), quad, true);
      const auto original(V);
      auto       b =
        assemble_lagrange_rhs(space, [](const Point &) { return 1.; }, 3);
      const auto                 integrals(b);
      libMesh::DenseVector<Real> s;
      V.lu_solve(b, s);
      Real charge = 0.;
      for (unsigned int i = 0; i < s.size(); ++i)
        {
          charge += s(i) * integrals(i);
          Real residual = -integrals(i);
          for (unsigned int j = 0; j < s.size(); ++j)
            residual += original(i, j) * s(j);
          REQUIRE(std::abs(residual) < 1e-11);
        }
      REQUIRE(charge / (4. * libMesh::pi) ==
              Catch::Approx(0.3667874).epsilon(0.02));
    }
}
