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
 * @file triangle_pair.cc
 * @brief Implementation of the neighboring type detection for triangles.
 */

#include "libmesh_bem/triangle_pair.h"

#include <libmesh/libmesh_common.h>

#include <vector>

namespace HierBEM
{
  namespace LibMeshBEM
  {
    using libMesh::Point;
    using libMesh::Real;

    PermutedTriangle::PermutedTriangle(
      const libMesh::Elem               &elem,
      const std::array<unsigned int, 3> &permutation)
      : perm(permutation)
    {
      libmesh_error_msg_if(elem.dim() != 2 || elem.n_vertices() != 3,
                           "Only triangular surface elements are supported, "
                           "but got element type "
                             << elem.type());

      for (unsigned int i = 0; i < 3; ++i)
        vertices[i] = elem.point(perm[i]);

      const Point n = (vertices[1] - vertices[0]).cross(vertices[2] - vertices[1]);
      jac           = n.norm();
      libmesh_error_msg_if(jac <= 0., "Degenerate triangle " << elem.id());

      // The normal is oriented with respect to the original vertex order. An
      // odd permutation flips the orientation.
      const Point n_orig =
        (elem.point(1) - elem.point(0)).cross(elem.point(2) - elem.point(0));
      unit_normal = n_orig / n_orig.norm();
    }


    Point
    PermutedTriangle::map_to_real(const Point &x_hat) const
    {
      return vertices[0] + x_hat(0) * (vertices[1] - vertices[0]) +
             x_hat(1) * (vertices[2] - vertices[1]);
    }


    Point
    PermutedTriangle::map_to_libmesh_reference(const Point &x_hat) const
    {
      // Barycentric coordinates with respect to the original vertex order.
      std::array<Real, 3> lambda;
      lambda[perm[0]] = 1. - x_hat(0);
      lambda[perm[1]] = x_hat(0) - x_hat(1);
      lambda[perm[2]] = x_hat(1);

      return Point(lambda[1], lambda[2], 0.);
    }


    TrianglePairInfo
    detect_triangle_pair(const libMesh::Elem &kx, const libMesh::Elem &ky)
    {
      TrianglePairInfo info;
      info.kx_permutation = {{0, 1, 2}};
      info.ky_permutation = {{0, 1, 2}};

      if (&kx == &ky || kx.id() == ky.id())
        {
          info.type = CellNeighboringType::SamePanel;
          return info;
        }

      // Collect shared vertices as pairs of local vertex indices.
      std::vector<std::pair<unsigned int, unsigned int>> shared;
      for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
          if (kx.node_id(i) == ky.node_id(j))
            shared.emplace_back(i, j);

      switch (shared.size())
        {
          case 0:
            info.type = CellNeighboringType::Regular;
            break;
          case 1:
            {
              const unsigned int i0 = shared[0].first;
              const unsigned int j0 = shared[0].second;
              info.type             = CellNeighboringType::CommonVertex;
              info.kx_permutation   = {{i0, (i0 + 1) % 3, (i0 + 2) % 3}};
              info.ky_permutation   = {{j0, (j0 + 1) % 3, (j0 + 2) % 3}};
              break;
            }
          case 2:
            {
              const auto [i0, j0]   = shared[0];
              const auto [i1, j1]   = shared[1];
              info.type             = CellNeighboringType::CommonEdge;
              info.kx_permutation   = {{i0, i1, 3 - i0 - i1}};
              info.ky_permutation   = {{j0, j1, 3 - j0 - j1}};
              break;
            }
          case 3:
            {
              // Two distinct elements with the same vertices: treat them as
              // the same panel with matching vertex order.
              info.type = CellNeighboringType::SamePanel;
              for (const auto &[i, j] : shared)
                info.ky_permutation[i] = j;
              break;
            }
          default:
            libmesh_error_msg("Invalid number of shared vertices");
        }

      return info;
    }
  } // namespace LibMeshBEM
} // namespace HierBEM
