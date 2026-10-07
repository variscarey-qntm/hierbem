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
 * @file galerkin_assembly.cc
 * @brief Implementation of the P0 Galerkin BEM assembly helpers.
 */

#include "libmesh_bem/galerkin_assembly.h"

namespace HierBEM
{
  namespace LibMeshBEM
  {
    using libMesh::Point;
    using libMesh::Real;

    std::vector<const libMesh::Elem *>
    collect_p0_elements(const libMesh::MeshBase &mesh)
    {
      std::vector<const libMesh::Elem *> elems;
      elems.reserve(mesh.n_active_elem());
      for (const libMesh::Elem *elem : mesh.active_element_ptr_range())
        {
          libmesh_error_msg_if(elem->dim() != 2 || elem->n_vertices() != 3,
                               "Only triangular surface elements are "
                               "supported, but got element type "
                                 << elem->type());
          elems.push_back(elem);
        }
      return elems;
    }


    libMesh::DenseVector<Real>
    assemble_p0_rhs(const std::vector<const libMesh::Elem *> &elems,
                    const std::function<Real(const Point &)> &g,
                    const unsigned int                        n_points)
    {
      std::vector<Point> qp;
      std::vector<Real>  qw;
      gauss_rule_on_sauter_reference_triangle(n_points, qp, qw);

      libMesh::DenseVector<Real> b(elems.size());
      for (std::size_t i = 0; i < elems.size(); ++i)
        {
          const PermutedTriangle t(*elems[i], {{0, 1, 2}});
          Real                   sum = 0.;
          for (std::size_t q = 0; q < qp.size(); ++q)
            sum += g(t.map_to_real(qp[q])) * qw[q];
          b(i) = sum * t.jacobian();
        }
      return b;
    }


    std::vector<Real>
    element_areas(const std::vector<const libMesh::Elem *> &elems)
    {
      std::vector<Real> areas(elems.size());
      for (std::size_t i = 0; i < elems.size(); ++i)
        areas[i] = PermutedTriangle(*elems[i], {{0, 1, 2}}).jacobian() / 2.;
      return areas;
    }
  } // namespace LibMeshBEM
} // namespace HierBEM
