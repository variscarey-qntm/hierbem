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
 * @file sauter_quadrature_triangle.h
 * @brief Sauter–Schwab quadrature for pairs of triangles, built on libMesh.
 *
 * This is the triangular counterpart of the quadrilateral Sauter quadrature
 * in @p quadrature/sauter_quadrature_tools.h, which is based on deal.II. The
 * four cell neighboring types (same panel, common edge, common vertex and
 * regular) are handled by the simplex coordinate transformations given in
 * S. A. Sauter and C. Schwab, "Boundary Element Methods", Springer 2011,
 * Section 5.2.1.
 *
 * All quadrature points are expressed in the Sauter reference triangle
 * \f$\hat{T} = \{(\hat{x}_1, \hat{x}_2) : 0 \le \hat{x}_2 \le \hat{x}_1 \le
 * 1\}\f$, with vertices \f$(0,0)\f$, \f$(1,0)\f$ and \f$(1,1)\f$. For a
 * triangle with (permuted) vertices \f$P_0, P_1, P_2\f$ the affine map is
 * \f$\chi(\hat{x}) = P_0 + \hat{x}_1 (P_1 - P_0) + \hat{x}_2 (P_2 - P_1)\f$.
 * In the singular cases, the shared vertex is mapped to \f$(0,0)\f$ and the
 * shared edge to \f$\hat{x}_2 = 0\f$ in both triangles.
 */

#ifndef HIERBEM_INCLUDE_LIBMESH_BEM_SAUTER_QUADRATURE_TRIANGLE_H_
#define HIERBEM_INCLUDE_LIBMESH_BEM_SAUTER_QUADRATURE_TRIANGLE_H_

#include <libmesh/point.h>

#include <array>
#include <vector>

namespace HierBEM
{
  namespace LibMeshBEM
  {
    /**
     * Different neighboring types of two triangles.
     */
    enum class CellNeighboringType
    {
      SamePanel,
      CommonEdge,
      CommonVertex,
      Regular
    };

    /**
     * Get the string representation of the cell neighboring type.
     */
    const char *
    cell_neighboring_type_name(const CellNeighboringType type);

    /**
     * Sauter quadrature order for the four cell neighboring types. As in the
     * deal.II based implementation, each order is the number of Gauss points
     * in each coordinate direction, i.e. a 1D rule of order @p n is exact for
     * polynomials up to degree \f$2n-1\f$. For the regular case, a libMesh
     * triangle Gauss rule of the same polynomial exactness \f$2n-1\f$ is used
     * on each triangle.
     */
    struct SauterQuadOrder
    {
      unsigned int same_panel    = 5;
      unsigned int common_edge   = 4;
      unsigned int common_vertex = 4;
      unsigned int regular       = 3;
    };

    /**
     * A quadrature point in the product space \f$\hat{T} \times \hat{T}\f$.
     * The weight already contains the Jacobian of the Sauter transformation.
     */
    struct QuadraturePointPair
    {
      libMesh::Point x_hat;
      libMesh::Point y_hat;
      libMesh::Real  weight;
    };

    /**
     * Sauter quadrature rules for triangle pairs. The 1D Gauss rules and the
     * triangle Gauss rule for the regular case are obtained from libMesh's
     * @p QGauss. All rules are precomputed in the constructor, so that the
     * sum of the weights of each rule is \f$|\hat{T}|^2 = 1/4\f$.
     */
    class SauterTriangleQuadrature
    {
    public:
      explicit SauterTriangleQuadrature(
        const SauterQuadOrder &quad_order = SauterQuadOrder());

      /**
       * Get the quadrature rule for the given cell neighboring type.
       */
      const std::vector<QuadraturePointPair> &
      get_rule(const CellNeighboringType type) const;

      const SauterQuadOrder &
      get_order() const
      {
        return quad_order;
      }

    private:
      SauterQuadOrder                  quad_order;
      std::vector<QuadraturePointPair> same_panel_rule;
      std::vector<QuadraturePointPair> common_edge_rule;
      std::vector<QuadraturePointPair> common_vertex_rule;
      std::vector<QuadraturePointPair> regular_rule;
    };

    /**
     * Gauss–Legendre rule with @p n_points points on \f$[0,1]\f$, obtained
     * from libMesh's 1D @p QGauss.
     */
    void
    gauss_rule_on_unit_interval(const unsigned int          n_points,
                                std::vector<libMesh::Real> &points,
                                std::vector<libMesh::Real> &weights);

    /**
     * Gauss rule on \f$\hat{T}\f$ with polynomial exactness \f$2n-1\f$,
     * obtained from libMesh's @p QGauss on the reference @p TRI3 and mapped
     * to \f$\hat{T}\f$.
     */
    void
    gauss_rule_on_sauter_reference_triangle(
      const unsigned int           n,
      std::vector<libMesh::Point> &points,
      std::vector<libMesh::Real>  &weights);

    /**
     * Number of subregions of the 4D unit cube for each neighboring type
     * (6, 5, 2 and 1 respectively).
     */
    unsigned int
    sauter_number_of_subregions(const CellNeighboringType type);

    /**
     * Transform the parametric coordinates \f$(\xi, \eta_1, \eta_2,
     * \eta_3)\f$ in the 4D unit cube to points in \f$\hat{T}\f$ for the
     * same panel case. @p k3_index is the subregion index in \f$[0, 6)\f$.
     *
     * @return Jacobian of the transformation.
     */
    libMesh::Real
    sauter_same_panel_parametric_coords_to_unit_cells(
      const std::array<libMesh::Real, 4> &parametric_coords,
      const unsigned int                  k3_index,
      libMesh::Point                     &kx_unit_cell_coords,
      libMesh::Point                     &ky_unit_cell_coords);

    /**
     * Same as above for the common edge case. @p k3_index is in \f$[0, 5)\f$.
     */
    libMesh::Real
    sauter_common_edge_parametric_coords_to_unit_cells(
      const std::array<libMesh::Real, 4> &parametric_coords,
      const unsigned int                  k3_index,
      libMesh::Point                     &kx_unit_cell_coords,
      libMesh::Point                     &ky_unit_cell_coords);

    /**
     * Same as above for the common vertex case. @p k3_index is in
     * \f$[0, 2)\f$.
     */
    libMesh::Real
    sauter_common_vertex_parametric_coords_to_unit_cells(
      const std::array<libMesh::Real, 4> &parametric_coords,
      const unsigned int                  k3_index,
      libMesh::Point                     &kx_unit_cell_coords,
      libMesh::Point                     &ky_unit_cell_coords);
  } // namespace LibMeshBEM
} // namespace HierBEM

#endif // HIERBEM_INCLUDE_LIBMESH_BEM_SAUTER_QUADRATURE_TRIANGLE_H_
