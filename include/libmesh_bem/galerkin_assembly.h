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
 * @file galerkin_assembly.h
 * @brief Galerkin BEM assembly on triangular libMesh surface meshes with
 * piecewise constant (P0) basis functions.
 */

#ifndef HIERBEM_INCLUDE_LIBMESH_BEM_GALERKIN_ASSEMBLY_H_
#define HIERBEM_INCLUDE_LIBMESH_BEM_GALERKIN_ASSEMBLY_H_

#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/elem.h>
#include <libmesh/mesh_base.h>

#include <functional>
#include <vector>

#include "libmesh_bem/sauter_quadrature_triangle.h"
#include "libmesh_bem/triangle_pair.h"

namespace HierBEM
{
  namespace LibMeshBEM
  {
    /**
     * Collect the active elements of @p mesh in a vector, which defines the
     * numbering of the P0 degrees of freedom: the i-th element in the
     * returned vector carries the i-th degree of freedom.
     */
    std::vector<const libMesh::Elem *>
    collect_p0_elements(const libMesh::MeshBase &mesh);

    /**
     * Assemble the dense Galerkin matrix
     * \f$A_{ij} = \int_{K_i} \int_{K_j} k(x, y) \, \mathrm{d}s_y \,
     * \mathrm{d}s_x\f$ for piecewise constant basis functions, where the
     * kernel is called as @p kernel(x, y, nx, ny).
     *
     * @param symmetric If true, only the upper triangle is computed and copied
     * to the lower triangle, which is valid for symmetric kernels such as the
     * single layer kernel.
     */
    template <typename Kernel>
    libMesh::DenseMatrix<libMesh::Real>
    assemble_p0_matrix(const std::vector<const libMesh::Elem *> &elems,
                       const Kernel                             &kernel,
                       const SauterTriangleQuadrature           &quad,
                       const bool                                symmetric)
    {
      const unsigned int                  n = elems.size();
      libMesh::DenseMatrix<libMesh::Real> A(n, n);

      for (unsigned int i = 0; i < n; ++i)
        for (unsigned int j = (symmetric ? i : 0); j < n; ++j)
          {
            const libMesh::Real a = integrate_on_triangle_pair(
              *elems[i],
              *elems[j],
              quad,
              [&kernel](const TriangleQuadraturePointData &x,
                        const TriangleQuadraturePointData &y) {
                return kernel(x.point, y.point, x.normal, y.normal);
              });

            A(i, j) = a;
            if (symmetric)
              A(j, i) = a;
          }

      return A;
    }

    /**
     * Assemble the Galerkin right hand side
     * \f$b_i = \int_{K_i} g(x) \, \mathrm{d}s_x\f$ for piecewise constant
     * test functions, using libMesh's triangle Gauss rule with @p n_points
     * points per direction (polynomial exactness \f$2n-1\f$).
     */
    libMesh::DenseVector<libMesh::Real>
    assemble_p0_rhs(
      const std::vector<const libMesh::Elem *>                   &elems,
      const std::function<libMesh::Real(const libMesh::Point &)> &g,
      const unsigned int                                          n_points);

    /**
     * Compute the areas of the triangles in @p elems.
     */
    std::vector<libMesh::Real>
    element_areas(const std::vector<const libMesh::Elem *> &elems);
  } // namespace LibMeshBEM
} // namespace HierBEM

#endif // HIERBEM_INCLUDE_LIBMESH_BEM_GALERKIN_ASSEMBLY_H_
