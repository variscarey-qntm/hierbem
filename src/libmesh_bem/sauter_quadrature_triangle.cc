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
 * @file sauter_quadrature_triangle.cc
 * @brief Implementation of the Sauter–Schwab quadrature for triangle pairs.
 */

#include "libmesh_bem/sauter_quadrature_triangle.h"

#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_order.h>
#include <libmesh/libmesh_common.h>
#include <libmesh/quadrature_gauss.h>

namespace HierBEM
{
  namespace LibMeshBEM
  {
    using libMesh::Point;
    using libMesh::Real;

    namespace
    {
      libMesh::Order
      exactness_order(const unsigned int n_points)
      {
        libmesh_error_msg_if(n_points == 0,
                             "The quadrature order must be positive.");
        return static_cast<libMesh::Order>(2 * n_points - 1);
      }

      /**
       * Build the 4D tensor product Gauss rule on the unit cube and push the
       * transformed points of all subregions into @p rule.
       */
      template <typename Transform>
      void
      build_singular_rule(const unsigned int                n_points,
                          const unsigned int                n_subregions,
                          Transform                       &&transform,
                          std::vector<QuadraturePointPair> &rule)
      {
        std::vector<Real> p, w;
        gauss_rule_on_unit_interval(n_points, p, w);

        const std::size_t n = p.size();
        rule.clear();
        rule.reserve(n * n * n * n * n_subregions);

        for (std::size_t i0 = 0; i0 < n; ++i0)
          for (std::size_t i1 = 0; i1 < n; ++i1)
            for (std::size_t i2 = 0; i2 < n; ++i2)
              for (std::size_t i3 = 0; i3 < n; ++i3)
                {
                  const std::array<Real, 4> coords{
                    {p[i0], p[i1], p[i2], p[i3]}};
                  const Real w4 = w[i0] * w[i1] * w[i2] * w[i3];

                  for (unsigned int k = 0; k < n_subregions; ++k)
                    {
                      QuadraturePointPair qp;
                      const Real jac = transform(coords, k, qp.x_hat, qp.y_hat);
                      qp.weight      = w4 * jac;
                      rule.push_back(qp);
                    }
                }
      }
    } // namespace


    const char *
    cell_neighboring_type_name(const CellNeighboringType type)
    {
      switch (type)
        {
          case CellNeighboringType::SamePanel:
            return "same panel";
          case CellNeighboringType::CommonEdge:
            return "common edge";
          case CellNeighboringType::CommonVertex:
            return "common vertex";
          case CellNeighboringType::Regular:
            return "regular";
        }
      return "unknown";
    }


    unsigned int
    sauter_number_of_subregions(const CellNeighboringType type)
    {
      switch (type)
        {
          case CellNeighboringType::SamePanel:
            return 6;
          case CellNeighboringType::CommonEdge:
            return 5;
          case CellNeighboringType::CommonVertex:
            return 2;
          case CellNeighboringType::Regular:
            return 1;
        }
      return 0;
    }


    void
    gauss_rule_on_unit_interval(const unsigned int n_points,
                                std::vector<Real> &points,
                                std::vector<Real> &weights)
    {
      libMesh::QGauss q(1, exactness_order(n_points));
      q.init(libMesh::EDGE2);

      points.resize(q.n_points());
      weights.resize(q.n_points());
      // Map from libMesh's reference interval \f$[-1,1]\f$ to \f$[0,1]\f$.
      for (unsigned int i = 0; i < q.n_points(); ++i)
        {
          points[i]  = 0.5 * (q.qp(i)(0) + 1.0);
          weights[i] = 0.5 * q.w(i);
        }
    }


    void
    gauss_rule_on_sauter_reference_triangle(const unsigned int  n,
                                            std::vector<Point> &points,
                                            std::vector<Real>  &weights)
    {
      libMesh::QGauss q(2, exactness_order(n));
      q.init(libMesh::TRI3);

      points.resize(q.n_points());
      weights.resize(q.n_points());
      // libMesh's reference triangle has the vertices (0,0), (1,0) and (0,1),
      // i.e. its reference coordinates are the barycentric coordinates
      // \f$(\lambda_1, \lambda_2)\f$. In \f$\hat{T}\f$, \f$\hat{x}_1 =
      // \lambda_1 + \lambda_2\f$ and \f$\hat{x}_2 = \lambda_2\f$. This map has
      // unit Jacobian, so that the weights are unchanged.
      for (unsigned int i = 0; i < q.n_points(); ++i)
        {
          points[i]  = Point(q.qp(i)(0) + q.qp(i)(1), q.qp(i)(1), 0.);
          weights[i] = q.w(i);
        }
    }


    Real
    sauter_same_panel_parametric_coords_to_unit_cells(
      const std::array<Real, 4> &parametric_coords,
      const unsigned int         k3_index,
      Point                     &kx_unit_cell_coords,
      Point                     &ky_unit_cell_coords)
    {
      const Real xi = parametric_coords[0];
      const Real e1 = parametric_coords[1];
      const Real e2 = parametric_coords[2];
      const Real e3 = parametric_coords[3];

      switch (k3_index)
        {
          case 0:
            kx_unit_cell_coords = Point(xi, xi * (1. - e1 + e1 * e2));
            ky_unit_cell_coords =
              Point(xi * (1. - e1 * e2 * e3), xi * (1. - e1));
            break;
          case 1:
            kx_unit_cell_coords =
              Point(xi * (1. - e1 * e2 * e3), xi * (1. - e1));
            ky_unit_cell_coords = Point(xi, xi * (1. - e1 + e1 * e2));
            break;
          case 2:
            kx_unit_cell_coords = Point(xi, xi * e1 * (1. - e2 + e2 * e3));
            ky_unit_cell_coords =
              Point(xi * (1. - e1 * e2), xi * e1 * (1. - e2));
            break;
          case 3:
            kx_unit_cell_coords =
              Point(xi * (1. - e1 * e2), xi * e1 * (1. - e2));
            ky_unit_cell_coords = Point(xi, xi * e1 * (1. - e2 + e2 * e3));
            break;
          case 4:
            kx_unit_cell_coords =
              Point(xi * (1. - e1 * e2 * e3), xi * e1 * (1. - e2 * e3));
            ky_unit_cell_coords = Point(xi, xi * e1 * (1. - e2));
            break;
          case 5:
            kx_unit_cell_coords = Point(xi, xi * e1 * (1. - e2));
            ky_unit_cell_coords =
              Point(xi * (1. - e1 * e2 * e3), xi * e1 * (1. - e2 * e3));
            break;
          default:
            libmesh_error_msg("Invalid subregion index " << k3_index
                                                         << " for same panel");
        }

      return xi * xi * xi * e1 * e1 * e2;
    }


    Real
    sauter_common_edge_parametric_coords_to_unit_cells(
      const std::array<Real, 4> &parametric_coords,
      const unsigned int         k3_index,
      Point                     &kx_unit_cell_coords,
      Point                     &ky_unit_cell_coords)
    {
      const Real xi = parametric_coords[0];
      const Real e1 = parametric_coords[1];
      const Real e2 = parametric_coords[2];
      const Real e3 = parametric_coords[3];

      switch (k3_index)
        {
          case 0:
            kx_unit_cell_coords = Point(xi, xi * e1 * e3);
            ky_unit_cell_coords =
              Point(xi * (1. - e1 * e2), xi * e1 * (1. - e2));
            return xi * xi * xi * e1 * e1;
          case 1:
            kx_unit_cell_coords = Point(xi, xi * e1);
            ky_unit_cell_coords =
              Point(xi * (1. - e1 * e2 * e3), xi * e1 * e2 * (1. - e3));
            break;
          case 2:
            kx_unit_cell_coords =
              Point(xi * (1. - e1 * e2), xi * e1 * (1. - e2));
            ky_unit_cell_coords = Point(xi, xi * e1 * e2 * e3);
            break;
          case 3:
            kx_unit_cell_coords =
              Point(xi * (1. - e1 * e2 * e3), xi * e1 * e2 * (1. - e3));
            ky_unit_cell_coords = Point(xi, xi * e1);
            break;
          case 4:
            kx_unit_cell_coords =
              Point(xi * (1. - e1 * e2 * e3), xi * e1 * (1. - e2 * e3));
            ky_unit_cell_coords = Point(xi, xi * e1 * e2);
            break;
          default:
            libmesh_error_msg("Invalid subregion index " << k3_index
                                                         << " for common edge");
        }

      return xi * xi * xi * e1 * e1 * e2;
    }


    Real
    sauter_common_vertex_parametric_coords_to_unit_cells(
      const std::array<Real, 4> &parametric_coords,
      const unsigned int         k3_index,
      Point                     &kx_unit_cell_coords,
      Point                     &ky_unit_cell_coords)
    {
      const Real xi = parametric_coords[0];
      const Real e1 = parametric_coords[1];
      const Real e2 = parametric_coords[2];
      const Real e3 = parametric_coords[3];

      switch (k3_index)
        {
          case 0:
            kx_unit_cell_coords = Point(xi, xi * e1);
            ky_unit_cell_coords = Point(xi * e2, xi * e2 * e3);
            break;
          case 1:
            kx_unit_cell_coords = Point(xi * e2, xi * e2 * e3);
            ky_unit_cell_coords = Point(xi, xi * e1);
            break;
          default:
            libmesh_error_msg("Invalid subregion index "
                              << k3_index << " for common vertex");
        }

      return xi * xi * xi * e2;
    }


    SauterTriangleQuadrature::SauterTriangleQuadrature(
      const SauterQuadOrder &quad_order_)
      : quad_order(quad_order_)
    {
      build_singular_rule(quad_order.same_panel,
                          6,
                          sauter_same_panel_parametric_coords_to_unit_cells,
                          same_panel_rule);
      build_singular_rule(quad_order.common_edge,
                          5,
                          sauter_common_edge_parametric_coords_to_unit_cells,
                          common_edge_rule);
      build_singular_rule(quad_order.common_vertex,
                          2,
                          sauter_common_vertex_parametric_coords_to_unit_cells,
                          common_vertex_rule);

      // The regular case is the tensor product of two triangle Gauss rules.
      std::vector<Point> p;
      std::vector<Real>  w;
      gauss_rule_on_sauter_reference_triangle(quad_order.regular, p, w);

      regular_rule.clear();
      regular_rule.reserve(p.size() * p.size());
      for (std::size_t i = 0; i < p.size(); ++i)
        for (std::size_t j = 0; j < p.size(); ++j)
          regular_rule.push_back(QuadraturePointPair{p[i], p[j], w[i] * w[j]});
    }


    const std::vector<QuadraturePointPair> &
    SauterTriangleQuadrature::get_rule(const CellNeighboringType type) const
    {
      switch (type)
        {
          case CellNeighboringType::SamePanel:
            return same_panel_rule;
          case CellNeighboringType::CommonEdge:
            return common_edge_rule;
          case CellNeighboringType::CommonVertex:
            return common_vertex_rule;
          case CellNeighboringType::Regular:
            return regular_rule;
        }
      libmesh_error_msg("Unknown cell neighboring type");
    }
  } // namespace LibMeshBEM
} // namespace HierBEM
