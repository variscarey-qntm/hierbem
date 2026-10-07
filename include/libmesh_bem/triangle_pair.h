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
 * @file triangle_pair.h
 * @brief Detection of the neighboring type of two libMesh triangles and the
 * corresponding affine maps from the Sauter reference triangle.
 */

#ifndef HIERBEM_INCLUDE_LIBMESH_BEM_TRIANGLE_PAIR_H_
#define HIERBEM_INCLUDE_LIBMESH_BEM_TRIANGLE_PAIR_H_

#include <libmesh/elem.h>
#include <libmesh/point.h>

#include <array>

#include "libmesh_bem/sauter_quadrature_triangle.h"

namespace HierBEM
{
  namespace LibMeshBEM
  {
    /**
     * Geometry of a flat triangle whose vertices are reordered by a local
     * vertex permutation, so that the shared vertex or edge of a singular
     * element pair is mapped to the same location in \f$\hat{T}\f$.
     */
    class PermutedTriangle
    {
    public:
      PermutedTriangle(const libMesh::Elem               &elem,
                       const std::array<unsigned int, 3> &permutation);

      /**
       * Map a point in the Sauter reference triangle \f$\hat{T}\f$ to real
       * space.
       */
      libMesh::Point
      map_to_real(const libMesh::Point &x_hat) const;

      /**
       * Map a point in the Sauter reference triangle \f$\hat{T}\f$ to the
       * reference coordinates of libMesh's @p TRI3 element (in the
       * <em>original</em> vertex order of the element), which can be passed
       * to libMesh's @p FEInterface::shape for evaluating shape functions.
       */
      libMesh::Point
      map_to_libmesh_reference(const libMesh::Point &x_hat) const;

      /**
       * Surface Jacobian of the map from \f$\hat{T}\f$ to the triangle, which
       * equals twice the triangle area.
       */
      libMesh::Real
      jacobian() const
      {
        return jac;
      }

      /**
       * Unit normal vector, oriented according to the original vertex order
       * of the libMesh element.
       */
      const libMesh::Point &
      normal() const
      {
        return unit_normal;
      }

    private:
      std::array<libMesh::Point, 3> vertices;
      std::array<unsigned int, 3>   perm;
      libMesh::Point                unit_normal;
      libMesh::Real                 jac;
    };

    /**
     * Result of the neighboring type detection for two triangles.
     */
    struct TrianglePairInfo
    {
      CellNeighboringType         type;
      std::array<unsigned int, 3> kx_permutation;
      std::array<unsigned int, 3> ky_permutation;
    };

    /**
     * Detect the neighboring type of two triangles in the same mesh by
     * comparing their vertex node ids and compute the vertex permutations
     * required by the Sauter quadrature:
     * - same panel: identical vertex order in both triangles;
     * - common edge: the shared edge is the first edge \f$P_0 P_1\f$ with the
     *   same orientation in both triangles;
     * - common vertex: the shared vertex is the first vertex \f$P_0\f$ in
     *   both triangles.
     *
     * Only the first three (vertex) nodes of the elements are considered.
     */
    TrianglePairInfo
    detect_triangle_pair(const libMesh::Elem &kx, const libMesh::Elem &ky);

    /**
     * Data at a quadrature point on one of the two triangles.
     */
    struct TriangleQuadraturePointData
    {
      /// Point in real space.
      libMesh::Point point;
      /// Point in libMesh's @p TRI3 reference coordinates.
      libMesh::Point reference_point;
      /// Unit normal vector.
      libMesh::Point normal;
    };

    /**
     * Loop over all Sauter quadrature point pairs for the two triangles
     * @p kx and @p ky and call
     * @p worker(x_data, y_data, JxW), where @p JxW is the quadrature weight
     * including the Jacobians of both triangles. This is the libMesh/triangle
     * counterpart of the cell-pair loops in the deal.II based implementation
     * and allows the evaluation of arbitrary shape functions via the
     * reference points.
     */
    template <typename Worker>
    void
    for_each_quadrature_point_pair(const libMesh::Elem            &kx,
                                   const libMesh::Elem            &ky,
                                   const SauterTriangleQuadrature &quad,
                                   Worker                        &&worker)
    {
      const TrianglePairInfo info = detect_triangle_pair(kx, ky);
      const PermutedTriangle tx(kx, info.kx_permutation);
      const PermutedTriangle ty(ky, info.ky_permutation);
      const libMesh::Real    jac = tx.jacobian() * ty.jacobian();

      TriangleQuadraturePointData x_data, y_data;
      x_data.normal = tx.normal();
      y_data.normal = ty.normal();

      for (const QuadraturePointPair &qp : quad.get_rule(info.type))
        {
          x_data.point           = tx.map_to_real(qp.x_hat);
          x_data.reference_point = tx.map_to_libmesh_reference(qp.x_hat);
          y_data.point           = ty.map_to_real(qp.y_hat);
          y_data.reference_point = ty.map_to_libmesh_reference(qp.y_hat);

          worker(x_data, y_data, qp.weight * jac);
        }
    }

    /**
     * Compute \f$\int_{K_x} \int_{K_y} f(x, y) \, \mathrm{d}s_y \,
     * \mathrm{d}s_x\f$ with the Sauter quadrature, where @p f is called as
     * @p f(x_data, y_data) and returns a scalar of type @p Number.
     */
    template <typename Number = libMesh::Real, typename Integrand>
    Number
    integrate_on_triangle_pair(const libMesh::Elem            &kx,
                               const libMesh::Elem            &ky,
                               const SauterTriangleQuadrature &quad,
                               Integrand                     &&f)
    {
      Number result(0);
      for_each_quadrature_point_pair(kx,
                                     ky,
                                     quad,
                                     [&](const TriangleQuadraturePointData &x,
                                         const TriangleQuadraturePointData &y,
                                         const libMesh::Real JxW) {
                                       result += f(x, y) * JxW;
                                     });
      return result;
    }
  } // namespace LibMeshBEM
} // namespace HierBEM

#endif // HIERBEM_INCLUDE_LIBMESH_BEM_TRIANGLE_PAIR_H_
