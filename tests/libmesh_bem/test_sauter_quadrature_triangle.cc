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
 * @file test_sauter_quadrature_triangle.cc
 * @brief Verify the Sauter–Schwab quadrature for triangle pairs.
 */

#include <libmesh/elem.h>
#include <libmesh/fe_map.h>
#include <libmesh/replicated_mesh.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <vector>

#include "libmesh_bem/galerkin_assembly.h"
#include "libmesh_bem/laplace_kernels.h"
#include "libmesh_bem/sauter_quadrature_triangle.h"
#include "libmesh_bem/triangle_pair.h"
#include "test_common.h"

using namespace HierBEM::LibMeshBEM;
using libMesh::Point;
using libMesh::Real;

namespace
{
  void
  add_triangle(libMesh::MeshBase                  &mesh,
               const std::array<unsigned int, 3> &nodes)
  {
    std::unique_ptr<libMesh::Elem> elem = libMesh::Elem::build(libMesh::TRI3);
    for (unsigned int i = 0; i < 3; ++i)
      elem->set_node(i, mesh.node_ptr(nodes[i]));
    mesh.add_elem(std::move(elem));
  }

  const std::array<Point, 3> vertices{
    {Point(0., 0., 0.), Point(1., 0.2, 0.), Point(0.3, 0.9, 0.1)}};

  /**
   * Mesh with a single triangle.
   */
  void
  build_single_triangle(libMesh::MeshBase &mesh)
  {
    for (unsigned int i = 0; i < 3; ++i)
      mesh.add_point(vertices[i], i);
    add_triangle(mesh, {{0, 1, 2}});
    mesh.prepare_for_use();
  }

  /**
   * Mesh with the four children of the uniform refinement of the triangle.
   */
  void
  build_refined_triangle(libMesh::MeshBase &mesh)
  {
    for (unsigned int i = 0; i < 3; ++i)
      mesh.add_point(vertices[i], i);
    mesh.add_point((vertices[0] + vertices[1]) / 2., 3);
    mesh.add_point((vertices[1] + vertices[2]) / 2., 4);
    mesh.add_point((vertices[2] + vertices[0]) / 2., 5);
    add_triangle(mesh, {{0, 3, 5}});
    add_triangle(mesh, {{3, 1, 4}});
    add_triangle(mesh, {{5, 4, 2}});
    add_triangle(mesh, {{4, 5, 3}});
    mesh.prepare_for_use();
  }

  Real
  slp_integral(const libMesh::Elem            &kx,
               const libMesh::Elem            &ky,
               const SauterTriangleQuadrature &quad)
  {
    const LaplaceSingleLayerKernel kernel;
    return integrate_on_triangle_pair(
      kx,
      ky,
      quad,
      [&kernel](const TriangleQuadraturePointData &x,
                const TriangleQuadraturePointData &y) {
        return kernel(x.point, y.point, x.normal, y.normal);
      });
  }
} // namespace


TEST_CASE("Sauter rules integrate constants exactly", "[libmesh_bem]")
{
  const SauterTriangleQuadrature quad;

  for (const CellNeighboringType type : {CellNeighboringType::SamePanel,
                                         CellNeighboringType::CommonEdge,
                                         CellNeighboringType::CommonVertex,
                                         CellNeighboringType::Regular})
    {
      INFO(cell_neighboring_type_name(type));
      Real sum = 0.;
      for (const QuadraturePointPair &qp : quad.get_rule(type))
        {
          sum += qp.weight;
          // All points must be inside the Sauter reference triangle.
          for (const Point *p : {&qp.x_hat, &qp.y_hat})
            {
              REQUIRE((*p)(1) >= -1e-14);
              REQUIRE((*p)(1) <= (*p)(0) + 1e-14);
              REQUIRE((*p)(0) <= 1. + 1e-14);
            }
        }
      REQUIRE(sum == Catch::Approx(0.25).epsilon(1e-13));
    }
}


TEST_CASE("Neighboring type detection and permutations", "[libmesh_bem]")
{
  libMesh::ReplicatedMesh mesh(
    HierBEM::LibMeshBEM::testing::libmesh_init->comm(), 3);
  build_refined_triangle(mesh);

  const libMesh::Elem &e0 = mesh.elem_ref(0);
  const libMesh::Elem &e1 = mesh.elem_ref(1);
  const libMesh::Elem &e3 = mesh.elem_ref(3);

  REQUIRE(detect_triangle_pair(e0, e0).type == CellNeighboringType::SamePanel);
  REQUIRE(detect_triangle_pair(e0, e1).type ==
          CellNeighboringType::CommonVertex);
  REQUIRE(detect_triangle_pair(e0, e3).type == CellNeighboringType::CommonEdge);

  // The shared vertex/edge must be mapped to the same physical points.
  for (const auto &[a, b] : {std::make_pair(&e0, &e1), std::make_pair(&e0, &e3)})
    {
      const TrianglePairInfo info = detect_triangle_pair(*a, *b);
      const PermutedTriangle tx(*a, info.kx_permutation);
      const PermutedTriangle ty(*b, info.ky_permutation);
      REQUIRE((tx.map_to_real(Point(0., 0.)) - ty.map_to_real(Point(0., 0.)))
                .norm() < 1e-14);
      if (info.type == CellNeighboringType::CommonEdge)
        REQUIRE(
          (tx.map_to_real(Point(1., 0.)) - ty.map_to_real(Point(1., 0.)))
            .norm() < 1e-14);
    }

  // The map to libMesh reference coordinates must be consistent with
  // libMesh's own geometric mapping.
  const std::array<std::array<unsigned int, 3>, 6> perms{{{{0, 1, 2}},
                                                          {{1, 2, 0}},
                                                          {{2, 0, 1}},
                                                          {{0, 2, 1}},
                                                          {{2, 1, 0}},
                                                          {{1, 0, 2}}}};
  for (const auto &perm : perms)
    {
      const PermutedTriangle t(e3, perm);
      REQUIRE(t.jacobian() == Catch::Approx(2. * e3.volume()));
      for (const Point x_hat : {Point(0.3, 0.1), Point(0.9, 0.7)})
        {
          const Point ref = t.map_to_libmesh_reference(x_hat);
          const Point p   = libMesh::FEMap::map(2, &e3, ref);
          REQUIRE((p - t.map_to_real(x_hat)).norm() < 1e-14);
        }
    }
}


TEST_CASE("Sauter quadrature is consistent under refinement", "[libmesh_bem]")
{
  // For the single layer kernel and a triangle T with children T_i of a
  // uniform refinement, \f$\int_T\int_T = \sum_{i,j} \int_{T_i}\int_{T_j}\f$.
  // This combines all four neighboring types. Moreover, by homogeneity of the
  // kernel, \f$\int_{T_i}\int_{T_i} = \frac{1}{8}\int_T\int_T\f$.
  libMesh::ReplicatedMesh coarse(
    HierBEM::LibMeshBEM::testing::libmesh_init->comm(), 3);
  libMesh::ReplicatedMesh fine(
    HierBEM::LibMeshBEM::testing::libmesh_init->comm(), 3);
  build_single_triangle(coarse);
  build_refined_triangle(fine);

  SauterQuadOrder order;
  order.same_panel    = 8;
  order.common_edge   = 8;
  order.common_vertex = 8;
  order.regular       = 8;
  const SauterTriangleQuadrature quad(order);

  const Real full = slp_integral(coarse.elem_ref(0), coarse.elem_ref(0), quad);

  Real sum = 0.;
  for (unsigned int i = 0; i < 4; ++i)
    {
      for (unsigned int j = 0; j < 4; ++j)
        sum += slp_integral(fine.elem_ref(i), fine.elem_ref(j), quad);

      REQUIRE(slp_integral(fine.elem_ref(i), fine.elem_ref(i), quad) ==
              Catch::Approx(full / 8.).epsilon(1e-12));
    }

  REQUIRE(sum == Catch::Approx(full).epsilon(1e-9));

  // The default (cheaper) orders must still be accurate.
  const SauterTriangleQuadrature default_quad;
  REQUIRE(slp_integral(coarse.elem_ref(0), coarse.elem_ref(0), default_quad) ==
          Catch::Approx(full).epsilon(1e-6));

  // Symmetry of the assembled single layer matrix.
  const auto elems = collect_p0_elements(fine);
  const auto V =
    assemble_p0_matrix(elems, LaplaceSingleLayerKernel(), quad, false);
  for (unsigned int i = 0; i < 4; ++i)
    for (unsigned int j = 0; j < 4; ++j)
      REQUIRE(V(i, j) == Catch::Approx(V(j, i)).epsilon(1e-10));
}


TEST_CASE("Analytical integral for the same panel case", "[libmesh_bem]")
{
  // Equilateral triangle T with unit edge length. The inner integral
  // \f$\int_T \frac{1}{|x-y|} \mathrm{d}s_y\f$ for \f$x \in T\f$ is
  // evaluated analytically: in polar coordinates around x it equals
  // \f$\sum_e d_e(x) \int_e \frac{1}{|x-y|} \mathrm{d}s_y\f$, where
  // \f$d_e(x)\f$ is the distance from x to the line of edge e, and the edge
  // integral is an @p asinh expression. The outer integral is computed with a
  // high-order triangle rule.
  libMesh::ReplicatedMesh mesh(
    HierBEM::LibMeshBEM::testing::libmesh_init->comm(), 3);
  mesh.add_point(Point(0., 0., 0.), 0);
  mesh.add_point(Point(1., 0., 0.), 1);
  mesh.add_point(Point(0.5, std::sqrt(3.) / 2., 0.), 2);
  add_triangle(mesh, {{0, 1, 2}});
  mesh.prepare_for_use();

  SauterQuadOrder order;
  order.same_panel = 10;
  const SauterTriangleQuadrature quad(order);
  const Real I = slp_integral(mesh.elem_ref(0), mesh.elem_ref(0), quad) * 4. *
                 libMesh::pi;

  const std::array<Point, 3> v{
    {mesh.point(0), mesh.point(1), mesh.point(2)}};
  const auto edge_potential = [&](const Point &x) {
    Real sum = 0.;
    for (unsigned int e = 0; e < 3; ++e)
      {
        const Point a = v[e];
        const Point b = v[(e + 1) % 3];
        const Point t = (b - a) / (b - a).norm();
        const Real  s0 = (a - x) * t;
        const Real  s1 = (b - x) * t;
        const Point foot = a - s0 * t;
        const Real  d    = (x - foot).norm();
        sum += d * (std::asinh(s1 / d) - std::asinh(s0 / d));
      }
    return sum;
  };

  std::vector<Point> qp;
  std::vector<Real>  qw;
  gauss_rule_on_sauter_reference_triangle(15, qp, qw);
  const PermutedTriangle t(mesh.elem_ref(0), {{0, 1, 2}});
  Real                   reference = 0.;
  for (std::size_t q = 0; q < qp.size(); ++q)
    reference += edge_potential(t.map_to_real(qp[q])) * qw[q] * t.jacobian();

  REQUIRE(I == Catch::Approx(reference).epsilon(1e-6));
}
