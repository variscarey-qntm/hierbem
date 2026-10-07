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
 * @file single_layer_entry.cc
 * @brief Implementation of single layer matrix entry evaluation.
 */

#include "libmesh_bem/single_layer_entry.h"

#include <libmesh/fe_interface.h>

#include "libmesh_bem/laplace_kernels.h"
#include "libmesh_bem/triangle_pair.h"

namespace HierBEM
{
  namespace LibMeshBEM
  {
    using libMesh::Point;
    using libMesh::Real;

    SingleLayerEntryEvaluator::SingleLayerEntryEvaluator(
      const std::vector<const libMesh::Elem *> &elems,
      const SauterTriangleQuadrature           &quad)
      : elems(elems)
      , supports(elems.size())
      , is_p0(true)
      , type(libMesh::CONSTANT, libMesh::MONOMIAL)
      , quad(&quad)
    {
      for (unsigned int i = 0; i < elems.size(); ++i)
        {
          libmesh_error_msg_if(!elems[i] || elems[i]->dim() != 2 ||
                                 elems[i]->n_vertices() != 3,
                               "Only triangular surface elements are "
                               "supported");
          supports[i].push_back({i, 0});
        }
    }


    SingleLayerEntryEvaluator::SingleLayerEntryEvaluator(
      const LagrangeTriangleSpace    &space,
      const SauterTriangleQuadrature &quad)
      : elems(space.elements())
      , supports(space.n_dofs())
      , is_p0(false)
      , type(space.fe_type())
      , quad(&quad)
    {
      for (unsigned int e = 0; e < elems.size(); ++e)
        {
          const auto &dofs = space.dof_indices(e);
          for (unsigned int a = 0; a < dofs.size(); ++a)
            supports[dofs[a]].push_back({e, a});
        }
    }


    Real
    SingleLayerEntryEvaluator::shape(const unsigned int elem,
                                     const unsigned int local_dof,
                                     const Point       &reference_point) const
    {
      if (is_p0)
        return 1.;
      return libMesh::FEInterface::shape(
        2, type, elems[elem], local_dof, reference_point);
    }


    Real
    SingleLayerEntryEvaluator::entry(const unsigned int i,
                                     const unsigned int j) const
    {
      libmesh_error_msg_if(i >= n_dofs() || j >= n_dofs(),
                           "DOF index (" << i << ", " << j
                                         << ") out of range for "
                                         << n_dofs() << " DOFs");
      const LaplaceSingleLayerKernel kernel;
      Real                           result = 0.;
      for (const Support &sx : supports[i])
        for (const Support &sy : supports[j])
          for_each_quadrature_point_pair(
            *elems[sx.elem],
            *elems[sy.elem],
            *quad,
            [&](const TriangleQuadraturePointData &x,
                const TriangleQuadraturePointData &y,
                const Real                         JxW) {
              result += kernel(x.point, y.point, x.normal, y.normal) * JxW *
                        shape(sx.elem, sx.local_dof, x.reference_point) *
                        shape(sy.elem, sy.local_dof, y.reference_point);
            });
      return result;
    }


    void
    butterflypack_single_layer_element(int *m, int *n, double *val, void *quant)
    {
      const auto &evaluator =
        *static_cast<const SingleLayerEntryEvaluator *>(quant);
      libmesh_error_msg_if(*m < 1 || *n < 1,
                           "ButterflyPACK indices must be one-based");
      *val = evaluator.entry(static_cast<unsigned int>(*m - 1),
                             static_cast<unsigned int>(*n - 1));
    }
  } // namespace LibMeshBEM
} // namespace HierBEM
