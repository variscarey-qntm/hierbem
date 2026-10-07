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
 * @brief Implementation of triangular Galerkin BEM assembly helpers.
 */

#include "libmesh_bem/galerkin_assembly.h"

#include <limits>

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


    LagrangeTriangleSpace::LagrangeTriangleSpace(const libMesh::MeshBase &mesh,
                                                 const libMesh::DofMap &dof_map)
      : elems(collect_p0_elements(mesh))
      , n(0)
    {
      libmesh_error_msg_if(!mesh.is_replicated(),
                           "Dense assembly requires a replicated mesh");
      libmesh_error_msg_if(dof_map.n_variables() != 1,
                           "A single Lagrange variable is required");
      type = dof_map.variable_type(0);
      libmesh_error_msg_if(type.family != libMesh::LAGRANGE ||
                             (type.order != libMesh::FIRST &&
                              type.order != libMesh::SECOND),
                           "Only FIRST and SECOND LAGRANGE are supported");
      libmesh_error_msg_if(dof_map.n_constrained_dofs() != 0,
                           "Constrained Lagrange spaces are not supported");
      libmesh_error_msg_if(dof_map.n_dofs() == 0 ||
                             dof_map.n_dofs() >
                               std::numeric_limits<unsigned int>::max(),
                           "Invalid dense system size");
      n = dof_map.n_dofs();
      dofs.resize(elems.size());
      for (unsigned int i = 0; i < elems.size(); ++i)
        {
          const auto &elem = *elems[i];
          libmesh_error_msg_if(
            (elem.type() != libMesh::TRI3 && elem.type() != libMesh::TRI6) ||
              (type.order == libMesh::SECOND && elem.type() != libMesh::TRI6) ||
              elem.p_level() != 0,
            "Use TRI3/TRI6 for P1 and TRI6 for P2, without p-refinement");
          if (elem.type() == libMesh::TRI6)
            for (unsigned int e = 0; e < 3; ++e)
              libmesh_error_msg_if(
                (elem.point(3 + e) -
                 (elem.point(e) + elem.point((e + 1) % 3)) / 2.)
                    .norm() >
                  1e-12 * (elem.point(e) - elem.point((e + 1) % 3)).norm(),
                "Only affine TRI6 triangles are supported");
          dof_map.dof_indices(&elem, dofs[i], 0);
          libmesh_error_msg_if(dofs[i].size() !=
                                 (type.order == libMesh::FIRST ? 3u : 6u),
                               "Unexpected number of Lagrange triangle DOFs");
          for (const auto d : dofs[i])
            libmesh_error_msg_if(d >= n, "Invalid Lagrange DOF index");
        }
    }


    libMesh::DenseVector<Real>
    assemble_lagrange_rhs(const LagrangeTriangleSpace              &space,
                          const std::function<Real(const Point &)> &g,
                          const unsigned int                        n_points)
    {
      std::vector<Point> qp;
      std::vector<Real>  qw;
      gauss_rule_on_sauter_reference_triangle(n_points, qp, qw);
      libMesh::DenseVector<Real> b(space.n_dofs());
      for (unsigned int i = 0; i < space.elements().size(); ++i)
        {
          const auto            &elem = *space.elements()[i];
          const auto            &dofs = space.dof_indices(i);
          const PermutedTriangle t(elem, {{0, 1, 2}});
          for (unsigned int q = 0; q < qp.size(); ++q)
            {
              const Point ref  = t.map_to_libmesh_reference(qp[q]);
              const Real value = g(t.map_to_real(qp[q])) * qw[q] * t.jacobian();
              for (unsigned int a = 0; a < dofs.size(); ++a)
                b(dofs[a]) += value * libMesh::FEInterface::shape(
                                        2, space.fe_type(), &elem, a, ref);
            }
        }
      return b;
    }
  } // namespace LibMeshBEM
} // namespace HierBEM
