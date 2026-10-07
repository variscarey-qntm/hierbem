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
 * @file single_layer_entry.h
 * @brief Evaluation of single entries of the Galerkin matrix of the Laplace
 * single layer potential for P0, P1 and P2 Lagrange elements on triangles,
 * e.g. for the element extraction interface of ButterflyPACK.
 */

#ifndef HIERBEM_INCLUDE_LIBMESH_BEM_SINGLE_LAYER_ENTRY_H_
#define HIERBEM_INCLUDE_LIBMESH_BEM_SINGLE_LAYER_ENTRY_H_

#include <libmesh/elem.h>
#include <libmesh/fe_type.h>

#include <vector>

#include "libmesh_bem/galerkin_assembly.h"
#include "libmesh_bem/sauter_quadrature_triangle.h"

namespace HierBEM
{
  namespace LibMeshBEM
  {
    /**
     * Computes a single entry
     * \f$V_{ij} = \int_\Gamma \int_\Gamma \frac{\phi_i(x) \phi_j(y)}
     * {4\pi|x-y|} \, \mathrm{d}s_y \, \mathrm{d}s_x\f$
     * of the assembled single layer matrix without assembling the matrix,
     * by summing the contributions of all element pairs in the supports of
     * \f$\phi_i\f$ and \f$\phi_j\f$. The result agrees with
     * assemble_p0_matrix() and assemble_lagrange_matrix() called with
     * <tt>symmetric = false</tt>.
     *
     * The DOF numbering is the one of collect_p0_elements() for P0, and the
     * libMesh DofMap numbering of LagrangeTriangleSpace for P1/P2.
     *
     * The mesh, the space and @p quad must outlive the evaluator. entry() is
     * const and uses no shared scratch data, so it may be called concurrently.
     */
    class SingleLayerEntryEvaluator
    {
    public:
      /**
       * Piecewise constant (P0) basis: DOF @p i is the indicator function of
       * @p elems[i].
       */
      SingleLayerEntryEvaluator(const std::vector<const libMesh::Elem *> &elems,
                                const SauterTriangleQuadrature           &quad);

      /**
       * Continuous FIRST (P1) or SECOND (P2) order Lagrange basis.
       */
      SingleLayerEntryEvaluator(const LagrangeTriangleSpace    &space,
                                const SauterTriangleQuadrature &quad);

      unsigned int
      n_dofs() const
      {
        return supports.size();
      }

      /**
       * Return the matrix entry \f$V_{ij}\f$ with zero-based DOF indices.
       */
      libMesh::Real
      entry(const unsigned int i, const unsigned int j) const;

    private:
      /// An element in the support of a basis function and the local index
      /// of the basis function on that element.
      struct Support
      {
        unsigned int elem;
        unsigned int local_dof;
      };

      libMesh::Real
      shape(const unsigned int    elem,
            const unsigned int    local_dof,
            const libMesh::Point &reference_point) const;

      std::vector<const libMesh::Elem *> elems;
      std::vector<std::vector<Support>>  supports;
      bool                               is_p0;
      libMesh::FEType                    type;
      const SauterTriangleQuadrature    *quad;
    };

    /**
     * Element extraction callback in the format of ButterflyPACK's
     * <tt>C_FuncZmn</tt>, as passed to
     * <tt>d_c_bpack_construct_element_compute</tt>: set <tt>*val</tt> to the
     * entry in row <tt>*m</tt> and column <tt>*n</tt>, which are one-based
     * indices in the original (natural, unpermuted) DOF ordering. @p quant is
     * the <tt>C2Fptr</tt> user pointer and must point to a
     * SingleLayerEntryEvaluator.
     */
    void
    butterflypack_single_layer_element(int    *m,
                                       int    *n,
                                       double *val,
                                       void   *quant);
  } // namespace LibMeshBEM
} // namespace HierBEM

#endif // HIERBEM_INCLUDE_LIBMESH_BEM_SINGLE_LAYER_ENTRY_H_
