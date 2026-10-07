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
 * @file solve-electrostatic-plate.cc
 * @brief Compute the surface charge density on a thin unit square plate held
 * at a constant electric potential, using a dense Galerkin BEM with P0,
 * P1 or P2 basis functions on a triangular mesh read with libMesh.
 *
 * The plate \f$\Gamma = [0,1]^2 \times \{0\}\f$ is held at the potential
 * \f$U\f$ in free space. The total surface charge density \f$\sigma\f$ (sum of
 * both sides of the plate) satisfies the first kind integral equation
 * \f[
 * \frac{1}{4\pi\varepsilon_0} \int_{\Gamma} \frac{\sigma(y)}{|x-y|}
 * \mathrm{d}s_y = U, \quad x \in \Gamma,
 * \f]
 * i.e. \f$V (\sigma / \varepsilon_0) = U\f$ with the Laplace single layer
 * operator \f$V\f$. The Galerkin matrix is assembled with the Sauter–Schwab
 * quadrature for triangles.
 *
 * Usage:
 * @code
 * solve-electrostatic-plate [--mesh unit_square.msh] [--potential 1.0]
 *   [--same-panel 5] [--common-edge 4] [--common-vertex 4] [--regular 3]
 *   [--order 1] [--output-prefix plate]
 * @endcode
 */

#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/libmesh.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/system.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>

#include "libmesh_bem/galerkin_assembly.h"
#include "libmesh_bem/gmsh_io.h"
#include "libmesh_bem/laplace_kernels.h"
#include "libmesh_bem/sauter_quadrature_triangle.h"

using namespace HierBEM::LibMeshBEM;
using libMesh::Real;

namespace
{
  // Vacuum permittivity in F/m.
  constexpr Real epsilon0 = 8.8541878128e-12;

  // Capacitance of the unit square plate in units of \f$4\pi\varepsilon_0
  // \cdot 1\,\mathrm{m}\f$, see e.g. C.-O. Hwang and M. Mascagni,
  // "Electrical capacitance of the unit cube", J. Appl. Phys. 95 (2004) and
  // references therein for the unit square.
  constexpr Real normalized_reference_capacitance = 0.3667874;

  /**
   * Write the mesh and the surface charge density in the Gmsh
   * MSH 2.2 ASCII format, which can be visualized directly in Gmsh.
   */
  void
  write_gmsh_element_data(const std::string                        &filename,
                          const libMesh::MeshBase                  &mesh,
                          const std::vector<const libMesh::Elem *> &elems,
                          const std::vector<std::vector<Real>>     &values,
                          const bool                                nodal)
  {
    std::ofstream out(filename);
    out << std::setprecision(16);
    out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";

    out << "$Nodes\n" << mesh.n_nodes() << "\n";
    for (const libMesh::Node *node : mesh.node_ptr_range())
      out << node->id() + 1 << " " << (*node)(0) << " " << (*node)(1) << " "
          << (*node)(2) << "\n";
    out << "$EndNodes\n";

    out << "$Elements\n" << elems.size() << "\n";
    for (std::size_t i = 0; i < elems.size(); ++i)
      {
        // Gmsh element types 2 and 9: 3-node and 6-node triangles.
        out << i + 1 << " " << (elems[i]->n_nodes() == 6 ? 9 : 2) << " 2 "
            << elems[i]->subdomain_id() << " " << elems[i]->subdomain_id();
        for (unsigned int v = 0; v < elems[i]->n_nodes(); ++v)
          out << " " << elems[i]->node_id(v) + 1;
        out << "\n";
      }
    out << "$EndElements\n";

    out << (nodal ? "$ElementNodeData" : "$ElementData")
        << "\n1\n\"surface charge density [C/m^2]\"\n1\n0.0\n3\n0\n1\n"
        << elems.size() << "\n";
    for (std::size_t i = 0; i < elems.size(); ++i)
      {
        out << i + 1;
        if (nodal)
          out << " " << values[i].size();
        for (const Real value : values[i])
          out << " " << value;
        out << "\n";
      }
    out << (nodal ? "$EndElementNodeData\n" : "$EndElementData\n");
  }
} // namespace


int
main(int argc, char **argv)
{
  libMesh::LibMeshInit init(argc, argv);

  const std::string mesh_file =
    libMesh::command_line_next("--mesh",
                               std::string(SOURCE_DIR "/unit_square.msh"));
  const Real         potential = libMesh::command_line_next("--potential", 1.0);
  const unsigned int degree    = libMesh::command_line_next("--order", 1u);
  libmesh_error_msg_if(degree > 2, "--order must be 0, 1 or 2");
  libmesh_error_msg_if(!std::isfinite(potential) || potential == 0.,
                       "--potential must be finite and nonzero");
  const std::string output_prefix =
    libMesh::command_line_next("--output-prefix", std::string("plate"));

  SauterQuadOrder order;
  order.same_panel =
    libMesh::command_line_next("--same-panel", order.same_panel);
  order.common_edge =
    libMesh::command_line_next("--common-edge", order.common_edge);
  order.common_vertex =
    libMesh::command_line_next("--common-vertex", order.common_vertex);
  order.regular = libMesh::command_line_next("--regular", order.regular);

  // Read the triangular mesh in the Gmsh MSH 4.1 binary format.
  libMesh::ReplicatedMesh mesh(init.comm(), 3);
  read_gmsh(mesh, mesh_file);
  if (degree == 2)
    mesh.all_second_order();
  mesh.print_info(libMesh::out);

  libMesh::EquationSystems               systems(mesh);
  std::unique_ptr<LagrangeTriangleSpace> space;
  if (degree != 0)
    {
      auto &system = systems.add_system<libMesh::System>("charge");
      system.add_variable("s",
                          degree == 1 ? libMesh::FIRST : libMesh::SECOND,
                          libMesh::LAGRANGE);
      systems.init();
      space =
        std::make_unique<LagrangeTriangleSpace>(mesh, system.get_dof_map());
    }

  const std::vector<const libMesh::Elem *> elems = collect_p0_elements(mesh);
  const std::vector<Real>                  areas = element_areas(elems);
  libmesh_error_msg_if(elems.empty(), "The surface mesh is empty");

  libMesh::out << "Sauter quadrature orders (same panel, common edge, "
                  "common vertex, regular): "
               << order.same_panel << ", " << order.common_edge << ", "
               << order.common_vertex << ", " << order.regular << std::endl;

  // Assemble the single layer matrix and the right hand side.
  const auto                     t0 = std::chrono::steady_clock::now();
  const SauterTriangleQuadrature quad(order);
  libMesh::DenseMatrix<Real>     V =
    space ?
          assemble_lagrange_matrix(*space, LaplaceSingleLayerKernel(), quad, true) :
          assemble_p0_matrix(elems, LaplaceSingleLayerKernel(), quad, true);
  const auto unit_potential = [](const libMesh::Point &) { return 1.; };
  const libMesh::DenseVector<Real> basis_integrals =
    space ? assemble_lagrange_rhs(*space, unit_potential, 3) :
            assemble_p0_rhs(elems, unit_potential, 2);
  libMesh::DenseVector<Real> b(basis_integrals);
  b.scale(potential);
  const auto t1 = std::chrono::steady_clock::now();

  // Solve \f$V s = b\f$ with \f$s = \sigma / \varepsilon_0\f$.
  libMesh::DenseVector<Real> s;
  V.lu_solve(b, s);
  const auto t2 = std::chrono::steady_clock::now();

  libMesh::DenseVector<Real> sigma(s);
  sigma.scale(epsilon0);

  Real total_charge = 0.;
  Real sigma_min    = sigma(0);
  Real sigma_max    = sigma(0);
  for (unsigned int i = 0; i < sigma.size(); ++i)
    {
      total_charge += sigma(i) * basis_integrals(i);
      sigma_min = std::min(sigma_min, sigma(i));
      sigma_max = std::max(sigma_max, sigma(i));
    }

  const Real capacitance = total_charge / potential;
  const Real reference_capacitance =
    normalized_reference_capacitance * 4. * libMesh::pi * epsilon0;

  libMesh::out << std::setprecision(8)
               << "Number of triangles: " << elems.size() << "\n"
               << "Basis order: P" << degree << "\n"
               << "Number of DOFs: " << s.size() << "\n"
               << "Assembly time [s]: "
               << std::chrono::duration<double>(t1 - t0).count() << "\n"
               << "Solve time [s]: "
               << std::chrono::duration<double>(t2 - t1).count() << "\n"
               << "Surface charge density range [C/m^2]: [" << sigma_min << ", "
               << sigma_max << "]\n"
               << "Total charge [C]: " << total_charge << "\n"
               << "Capacitance [F]: " << capacitance << "\n"
               << "Normalized capacitance C/(4 pi eps0): "
               << capacitance / (4. * libMesh::pi * epsilon0) << "\n"
               << "Reference capacitance [F]: " << reference_capacitance << "\n"
               << "Relative error: "
               << std::abs(capacitance - reference_capacitance) /
                    reference_capacitance
               << std::endl;

  if (init.comm().rank() == 0)
    {
      std::ofstream csv(output_prefix + "-charge-density.csv");
      csv << std::setprecision(16) << "x,y,z,area,sigma\n";
      std::vector<std::vector<Real>>      values(elems.size());
      const std::array<libMesh::Point, 6> reference_nodes{
        {libMesh::Point(0., 0.),
         libMesh::Point(1., 0.),
         libMesh::Point(0., 1.),
         libMesh::Point(0.5, 0.),
         libMesh::Point(0.5, 0.5),
         libMesh::Point(0., 0.5)}};
      for (std::size_t i = 0; i < elems.size(); ++i)
        {
          const auto evaluate = [&](const libMesh::Point &ref) {
            if (!space)
              return sigma(i);
            Real        value = 0.;
            const auto &dofs  = space->dof_indices(i);
            for (unsigned int a = 0; a < dofs.size(); ++a)
              value +=
                sigma(dofs[a]) * libMesh::FEInterface::shape(
                                   2, space->fe_type(), elems[i], a, ref);
            return value;
          };
          const libMesh::Point c = elems[i]->vertex_average();
          csv << c(0) << "," << c(1) << "," << c(2) << "," << areas[i] << ","
              << evaluate(libMesh::Point(1. / 3., 1. / 3.)) << "\n";
          if (space)
            for (unsigned int a = 0; a < elems[i]->n_nodes(); ++a)
              values[i].push_back(evaluate(reference_nodes[a]));
          else
            values[i].push_back(sigma(i));
        }

      write_gmsh_element_data(output_prefix + "-charge-density.msh",
                              mesh,
                              elems,
                              values,
                              degree != 0);
      libMesh::out << "Results written to " << output_prefix
                   << "-charge-density.{csv,msh}" << std::endl;
    }

  return 0;
}
