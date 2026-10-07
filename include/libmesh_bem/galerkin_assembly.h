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
 * piecewise constant (P0) and continuous Lagrange (P1/P2) basis functions.
 */

#ifndef HIERBEM_INCLUDE_LIBMESH_BEM_GALERKIN_ASSEMBLY_H_
#define HIERBEM_INCLUDE_LIBMESH_BEM_GALERKIN_ASSEMBLY_H_

#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/fe_interface.h>
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

    /**
     * A single-variable FIRST or SECOND LAGRANGE space on replicated flat
     * triangle meshes.
     * The libMesh DofMap must be initialized, without hanging-node constraints.
     * SECOND order requires TRI6 elements with nodes at edge midpoints.
     * The mesh must outlive the space and must not change after construction.
     */
    class LagrangeTriangleSpace
    {
    public:
      LagrangeTriangleSpace(const libMesh::MeshBase &mesh,
                            const libMesh::DofMap   &dof_map);

      const std::vector<const libMesh::Elem *> &
      elements() const
      {
        return elems;
      }

      const std::vector<libMesh::dof_id_type> &
      dof_indices(const unsigned int i) const
      {
        return dofs[i];
      }

      const libMesh::FEType &
      fe_type() const
      {
        return type;
      }

      unsigned int
      n_dofs() const
      {
        return n;
      }

    private:
      std::vector<const libMesh::Elem *>             elems;
      std::vector<std::vector<libMesh::dof_id_type>> dofs;
      libMesh::FEType                                type;
      unsigned int                                   n;
    };

    /**
     * Assemble the dense matrix with libMesh Lagrange trial and test functions.
     * Symmetric kernels may reuse transposed off-diagonal element blocks.
     */
    template <typename Kernel>
    libMesh::DenseMatrix<libMesh::Real>
    assemble_lagrange_matrix(const LagrangeTriangleSpace    &space,
                             const Kernel                   &kernel,
                             const SauterTriangleQuadrature &quad,
                             const bool                      symmetric)
    {
      libMesh::DenseMatrix<libMesh::Real> A(space.n_dofs(), space.n_dofs());
      const auto                         &elems = space.elements();
      for (unsigned int i = 0; i < elems.size(); ++i)
        for (unsigned int j = (symmetric ? i : 0); j < elems.size(); ++j)
          {
            const auto                         &dx = space.dof_indices(i);
            const auto                         &dy = space.dof_indices(j);
            libMesh::DenseMatrix<libMesh::Real> block(dx.size(), dy.size());
            std::vector<libMesh::Real>          px(dx.size()), py(dy.size());
            for_each_quadrature_point_pair(
              *elems[i],
              *elems[j],
              quad,
              [&](const TriangleQuadraturePointData &x,
                  const TriangleQuadraturePointData &y,
                  const libMesh::Real                JxW) {
                for (unsigned int a = 0; a < dx.size(); ++a)
                  px[a] = libMesh::FEInterface::shape(
                    2, space.fe_type(), elems[i], a, x.reference_point);
                for (unsigned int b = 0; b < dy.size(); ++b)
                  py[b] = libMesh::FEInterface::shape(
                    2, space.fe_type(), elems[j], b, y.reference_point);
                const libMesh::Real value =
                  kernel(x.point, y.point, x.normal, y.normal) * JxW;
                for (unsigned int a = 0; a < dx.size(); ++a)
                  for (unsigned int b = 0; b < dy.size(); ++b)
                    block(a, b) += value * px[a] * py[b];
              });
            for (unsigned int a = 0; a < dx.size(); ++a)
              for (unsigned int b = 0; b < dy.size(); ++b)
                {
                  A(dx[a], dy[b]) += block(a, b);
                  if (symmetric && i != j)
                    A(dy[b], dx[a]) += block(a, b);
                }
          }
      return A;
    }

    /**
     * Assemble b_i = integral g(x) phi_i(x) ds using the triangle Gauss rule.
     * With g=1, the result also supplies basis integrals for total charge.
     */
    libMesh::DenseVector<libMesh::Real>
    assemble_lagrange_rhs(
      const LagrangeTriangleSpace                                &space,
      const std::function<libMesh::Real(const libMesh::Point &)> &g,
      const unsigned int                                          n_points);
  } // namespace LibMeshBEM
} // namespace HierBEM

#endif // HIERBEM_INCLUDE_LIBMESH_BEM_GALERKIN_ASSEMBLY_H_
