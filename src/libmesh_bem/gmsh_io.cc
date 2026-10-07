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
 * @file gmsh_io.cc
 * @brief Implementation of the Gmsh MSH 4.1 binary reader.
 */

#include "libmesh_bem/gmsh_io.h"

#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/gmsh_io.h>
#include <libmesh/libmesh_common.h>
#include <libmesh/node.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace HierBEM
{
  namespace LibMeshBEM
  {
    namespace
    {
      /**
       * Low level binary reader handling the byte order and the size of
       * @p size_t values recorded in the file header.
       */
      class BinaryReader
      {
      public:
        BinaryReader(std::istream &in_, const std::string &filename_)
          : in(in_)
          , filename(filename_)
        {}

        void
        set_swap(const bool swap)
        {
          swap_bytes = swap;
        }

        void
        set_size_t_bytes(const unsigned int n)
        {
          libmesh_error_msg_if(n != 4 && n != 8,
                               "Unsupported data size " << n << " in "
                                                        << filename);
          size_t_bytes = n;
        }

        std::int32_t
        read_int()
        {
          std::int32_t v;
          read_raw(&v, sizeof(v));
          return v;
        }

        double
        read_double()
        {
          double v;
          read_raw(&v, sizeof(v));
          return v;
        }

        std::uint64_t
        read_size_t()
        {
          if (size_t_bytes == 8)
            {
              std::uint64_t v;
              read_raw(&v, sizeof(v));
              return v;
            }
          std::uint32_t v;
          read_raw(&v, sizeof(v));
          return v;
        }

        void
        skip(const std::uint64_t n_bytes)
        {
          in.ignore(static_cast<std::streamsize>(n_bytes));
          check();
        }

        unsigned int
        size_t_size() const
        {
          return size_t_bytes;
        }

      private:
        void
        read_raw(void *dst, const std::size_t n)
        {
          in.read(static_cast<char *>(dst), static_cast<std::streamsize>(n));
          check();
          if (swap_bytes)
            {
              char *c = static_cast<char *>(dst);
              std::reverse(c, c + n);
            }
        }

        void
        check() const
        {
          libmesh_error_msg_if(!in,
                               "Unexpected end of file or read error in "
                                 << filename);
        }

        std::istream     &in;
        const std::string filename;
        bool              swap_bytes   = false;
        unsigned int      size_t_bytes = 8;
      };


      /**
       * Number of nodes of Gmsh element types, used for skipping elements
       * which are not imported.
       */
      unsigned int
      gmsh_element_n_nodes(const int type)
      {
        static const std::map<int, unsigned int> n_nodes = {
          {1, 2},   {2, 3},   {3, 4},   {4, 4},   {5, 8},   {6, 6},
          {7, 5},   {8, 3},   {9, 6},   {10, 9},  {11, 10}, {12, 27},
          {13, 18}, {14, 14}, {15, 1},  {16, 8},  {17, 20}, {18, 15},
          {19, 13}, {20, 9},  {21, 10}, {22, 12}, {23, 15}, {24, 15},
          {25, 21}, {26, 4},  {27, 5},  {28, 6},  {29, 20}, {30, 35},
          {31, 56}, {36, 16}, {37, 25}, {38, 36}, {92, 64}, {93, 125}};

        const auto it = n_nodes.find(type);
        libmesh_error_msg_if(it == n_nodes.end(),
                             "Unsupported Gmsh element type " << type);
        return it->second;
      }


      /**
       * Map supported Gmsh 2D element types to libMesh element types.
       */
      bool
      gmsh_to_libmesh_type(const int type, libMesh::ElemType &libmesh_type)
      {
        switch (type)
          {
            case 2:
              libmesh_type = libMesh::TRI3;
              return true;
            case 3:
              libmesh_type = libMesh::QUAD4;
              return true;
            case 9:
              libmesh_type = libMesh::TRI6;
              return true;
            default:
              return false;
          }
      }


      /**
       * Skip lines until the line @p end_tag is found.
       */
      void
      skip_to(std::istream      &in,
              const std::string &end_tag,
              const std::string &filename)
      {
        std::string line;
        while (std::getline(in, line))
          {
            if (!line.empty() && line.back() == '\r')
              line.pop_back();
            if (line == end_tag)
              return;
          }
        libmesh_error_msg("Cannot find " << end_tag << " in " << filename);
      }


      struct MeshFormat
      {
        double version   = 0.;
        int    file_type = 0;
        int    data_size = 0;
      };


      /**
       * Read the @p $MeshFormat section; the stream must be positioned
       * directly after the @p $MeshFormat line.
       */
      MeshFormat
      read_mesh_format(std::istream &in, const std::string &filename)
      {
        MeshFormat  format;
        std::string line;
        std::getline(in, line);
        std::istringstream ss(line);
        ss >> format.version >> format.file_type >> format.data_size;
        libmesh_error_msg_if(!ss, "Invalid $MeshFormat in " << filename);
        return format;
      }


      std::string
      trimmed(std::string s)
      {
        while (!s.empty() &&
               (s.back() == '\r' || s.back() == ' ' || s.back() == '\t'))
          s.pop_back();
        return s;
      }
    } // namespace


    void
    read_gmsh41_binary(libMesh::MeshBase &mesh, const std::string &filename)
    {
      std::ifstream in(filename, std::ios::in | std::ios::binary);
      libmesh_error_msg_if(!in, "Cannot open the mesh file " << filename);

      BinaryReader reader(in, filename);

      // Entity tag of surfaces => physical tag.
      std::map<int, int>                                surface_phys;
      std::map<int, std::string>                        phys_names;
      std::unordered_map<std::uint64_t, libMesh::Point> node_coords;
      bool                                              has_format = false;
      bool                                              has_nodes  = false;
      bool                                              has_elems  = false;

      struct ElementRecord
      {
        libMesh::ElemType          type;
        int                        subdomain;
        std::vector<std::uint64_t> nodes;
      };
      std::vector<ElementRecord> elements;

      std::string line;
      while (std::getline(in, line))
        {
          line = trimmed(line);
          if (line.empty())
            continue;

          if (line == "$MeshFormat")
            {
              const MeshFormat format = read_mesh_format(in, filename);
              libmesh_error_msg_if(format.version < 4.1 || format.version >= 5.,
                                   "Only Gmsh MSH 4.1 is supported, but "
                                     << filename << " has version "
                                     << format.version);
              libmesh_error_msg_if(format.file_type != 1,
                                   filename << " is not a binary MSH file");
              reader.set_size_t_bytes(format.data_size);

              // Endianness detection: the integer 1 is written in binary.
              std::int32_t one;
              in.read(reinterpret_cast<char *>(&one), sizeof(one));
              if (one != 1)
                {
                  char *c = reinterpret_cast<char *>(&one);
                  std::reverse(c, c + sizeof(one));
                  libmesh_error_msg_if(one != 1,
                                       "Invalid endianness marker in "
                                         << filename);
                  reader.set_swap(true);
                }
              skip_to(in, "$EndMeshFormat", filename);
              has_format = true;
            }
          else if (line == "$PhysicalNames")
            {
              // This section is always written in ASCII.
              std::getline(in, line);
              const int n = std::stoi(line);
              for (int i = 0; i < n; ++i)
                {
                  std::getline(in, line);
                  std::istringstream ss(line);
                  int                dim, tag;
                  ss >> dim >> tag;
                  std::string name;
                  std::getline(ss, name);
                  const auto first = name.find('"');
                  const auto last  = name.rfind('"');
                  if (first != std::string::npos && last > first)
                    name = name.substr(first + 1, last - first - 1);
                  if (dim == 2)
                    phys_names[tag] = name;
                }
              skip_to(in, "$EndPhysicalNames", filename);
            }
          else if (line == "$Entities")
            {
              libmesh_error_msg_if(!has_format,
                                   "Missing $MeshFormat in " << filename);
              const std::uint64_t n_points   = reader.read_size_t();
              const std::uint64_t n_curves   = reader.read_size_t();
              const std::uint64_t n_surfaces = reader.read_size_t();
              const std::uint64_t n_volumes  = reader.read_size_t();

              for (std::uint64_t i = 0; i < n_points; ++i)
                {
                  reader.read_int(); // tag
                  reader.skip(3 * sizeof(double));
                  const std::uint64_t n_phys = reader.read_size_t();
                  reader.skip(n_phys * sizeof(std::int32_t));
                }

              const auto read_entities = [&](const std::uint64_t n,
                                             const bool          is_surface) {
                for (std::uint64_t i = 0; i < n; ++i)
                  {
                    const int tag = reader.read_int();
                    reader.skip(6 * sizeof(double)); // bounding box
                    const std::uint64_t n_phys = reader.read_size_t();
                    for (std::uint64_t p = 0; p < n_phys; ++p)
                      {
                        const int phys = reader.read_int();
                        // Only the first physical tag is used.
                        if (is_surface && p == 0)
                          surface_phys[tag] = std::abs(phys);
                      }
                    const std::uint64_t n_bounding = reader.read_size_t();
                    reader.skip(n_bounding * sizeof(std::int32_t));
                  }
              };
              read_entities(n_curves, false);
              read_entities(n_surfaces, true);
              read_entities(n_volumes, false);

              skip_to(in, "$EndEntities", filename);
            }
          else if (line == "$PartitionedEntities")
            {
              libmesh_error_msg(
                "Partitioned Gmsh meshes are not supported: " << filename);
            }
          else if (line == "$Nodes")
            {
              libmesh_error_msg_if(!has_format,
                                   "Missing $MeshFormat in " << filename);
              const std::uint64_t n_blocks = reader.read_size_t();
              const std::uint64_t n_nodes  = reader.read_size_t();
              reader.read_size_t(); // min node tag
              reader.read_size_t(); // max node tag
              node_coords.reserve(n_nodes);

              std::vector<std::uint64_t> tags;
              for (std::uint64_t b = 0; b < n_blocks; ++b)
                {
                  const int entity_dim = reader.read_int();
                  reader.read_int(); // entity tag
                  const int           parametric = reader.read_int();
                  const std::uint64_t n_in_block = reader.read_size_t();

                  tags.resize(n_in_block);
                  for (auto &t : tags)
                    t = reader.read_size_t();

                  const unsigned int n_param =
                    parametric ? static_cast<unsigned int>(entity_dim) : 0;
                  for (std::uint64_t i = 0; i < n_in_block; ++i)
                    {
                      const double x = reader.read_double();
                      const double y = reader.read_double();
                      const double z = reader.read_double();
                      reader.skip(n_param * sizeof(double));
                      node_coords[tags[i]] = libMesh::Point(x, y, z);
                    }
                }
              skip_to(in, "$EndNodes", filename);
              has_nodes = true;
            }
          else if (line == "$Elements")
            {
              libmesh_error_msg_if(!has_format,
                                   "Missing $MeshFormat in " << filename);
              const std::uint64_t n_blocks = reader.read_size_t();
              reader.read_size_t(); // number of elements
              reader.read_size_t(); // min element tag
              reader.read_size_t(); // max element tag

              for (std::uint64_t b = 0; b < n_blocks; ++b)
                {
                  const int           entity_dim = reader.read_int();
                  const int           entity_tag = reader.read_int();
                  const int           gmsh_type  = reader.read_int();
                  const std::uint64_t n_in_block = reader.read_size_t();
                  const unsigned int  n_vertices =
                    gmsh_element_n_nodes(gmsh_type);

                  libMesh::ElemType libmesh_type;
                  const bool        import =
                    entity_dim == 2 &&
                    gmsh_to_libmesh_type(gmsh_type, libmesh_type);
                  libmesh_error_msg_if(entity_dim == 2 && !import,
                                       "Unsupported 2D Gmsh element type "
                                         << gmsh_type << " in " << filename);

                  if (!import)
                    {
                      reader.skip(n_in_block * (1 + n_vertices) *
                                  reader.size_t_size());
                      continue;
                    }

                  const auto it_phys   = surface_phys.find(entity_tag);
                  const int  subdomain = it_phys != surface_phys.end() ?
                                           it_phys->second :
                                           entity_tag;

                  for (std::uint64_t e = 0; e < n_in_block; ++e)
                    {
                      reader.read_size_t(); // element tag
                      ElementRecord rec;
                      rec.type      = libmesh_type;
                      rec.subdomain = subdomain;
                      rec.nodes.resize(n_vertices);
                      for (auto &n : rec.nodes)
                        n = reader.read_size_t();
                      elements.push_back(std::move(rec));
                    }
                }
              skip_to(in, "$EndElements", filename);
              has_elems = true;
              // All required data are available.
              break;
            }
          else if (line.size() > 1 && line[0] == '$' &&
                   line.compare(0, 4, "$End") != 0)
            {
              // Skip unknown or unused sections.
              skip_to(in, "$End" + line.substr(1), filename);
            }
        }

      libmesh_error_msg_if(!has_nodes || !has_elems,
                           "Missing $Nodes or $Elements section in "
                             << filename);
      libmesh_error_msg_if(elements.empty(),
                           "No supported 2D elements found in " << filename);

      // Build the libMesh mesh. Only nodes referenced by imported elements
      // are added.
      mesh.clear();
      mesh.set_mesh_dimension(2);

      std::unordered_map<std::uint64_t, libMesh::Node *> tag_to_node;
      libMesh::dof_id_type                               next_node_id = 0;
      libMesh::dof_id_type                               next_elem_id = 0;

      for (const ElementRecord &rec : elements)
        {
          std::unique_ptr<libMesh::Elem> elem = libMesh::Elem::build(rec.type);
          for (unsigned int i = 0; i < rec.nodes.size(); ++i)
            {
              auto it = tag_to_node.find(rec.nodes[i]);
              if (it == tag_to_node.end())
                {
                  const auto c = node_coords.find(rec.nodes[i]);
                  libmesh_error_msg_if(c == node_coords.end(),
                                       "Unknown node tag "
                                         << rec.nodes[i] << " in " << filename);
                  libMesh::Node *node =
                    mesh.add_point(c->second, next_node_id++);
                  it = tag_to_node.emplace(rec.nodes[i], node).first;
                }
              elem->set_node(i) = it->second;
            }
          elem->subdomain_id() =
            static_cast<libMesh::subdomain_id_type>(rec.subdomain);
          elem->set_id(next_elem_id++);
          mesh.add_elem(std::move(elem));
        }

      for (const auto &[tag, name] : phys_names)
        mesh.subdomain_name(static_cast<libMesh::subdomain_id_type>(tag)) =
          name;

      mesh.prepare_for_use();
    }


    void
    read_gmsh(libMesh::MeshBase &mesh, const std::string &filename)
    {
      std::ifstream in(filename, std::ios::in | std::ios::binary);
      libmesh_error_msg_if(!in, "Cannot open the mesh file " << filename);

      std::string line;
      while (std::getline(in, line))
        if (trimmed(line) == "$MeshFormat")
          {
            const MeshFormat format = read_mesh_format(in, filename);
            in.close();
            if (format.file_type == 1)
              read_gmsh41_binary(mesh, filename);
            else
              {
                libMesh::GmshIO gmsh_io(mesh);
                gmsh_io.read(filename);
                mesh.prepare_for_use();
              }
            return;
          }

      libmesh_error_msg("Missing $MeshFormat in " << filename);
    }
  } // namespace LibMeshBEM
} // namespace HierBEM
