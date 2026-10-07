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
 * @file laplace_kernels.h
 * @brief Kernel functions of the Laplace boundary integral operators in 3D,
 * using libMesh point types.
 */

#ifndef HIERBEM_INCLUDE_LIBMESH_BEM_LAPLACE_KERNELS_H_
#define HIERBEM_INCLUDE_LIBMESH_BEM_LAPLACE_KERNELS_H_

#include <libmesh/libmesh_common.h>
#include <libmesh/point.h>

namespace HierBEM
{
  namespace LibMeshBEM
  {
    /**
     * Single layer kernel \f$\frac{1}{4\pi|x-y|}\f$.
     */
    struct LaplaceSingleLayerKernel
    {
      libMesh::Real
      operator()(const libMesh::Point &x,
                 const libMesh::Point &y,
                 const libMesh::Point & /* nx */,
                 const libMesh::Point & /* ny */) const
      {
        return 1.0 / (4.0 * libMesh::pi * (x - y).norm());
      }
    };

    /**
     * Double layer kernel \f$\frac{\langle x-y, n_y \rangle}{4\pi|x-y|^3}\f$.
     */
    struct LaplaceDoubleLayerKernel
    {
      libMesh::Real
      operator()(const libMesh::Point &x,
                 const libMesh::Point &y,
                 const libMesh::Point & /* nx */,
                 const libMesh::Point &ny) const
      {
        const libMesh::Point d = x - y;
        const libMesh::Real  r = d.norm();
        return (d * ny) / (4.0 * libMesh::pi * r * r * r);
      }
    };

    /**
     * Adjoint double layer kernel
     * \f$\frac{\langle y-x, n_x \rangle}{4\pi|x-y|^3}\f$.
     */
    struct LaplaceAdjointDoubleLayerKernel
    {
      libMesh::Real
      operator()(const libMesh::Point &x,
                 const libMesh::Point &y,
                 const libMesh::Point &nx,
                 const libMesh::Point & /* ny */) const
      {
        const libMesh::Point d = y - x;
        const libMesh::Real  r = d.norm();
        return (d * nx) / (4.0 * libMesh::pi * r * r * r);
      }
    };
  } // namespace LibMeshBEM
} // namespace HierBEM

#endif // HIERBEM_INCLUDE_LIBMESH_BEM_LAPLACE_KERNELS_H_
