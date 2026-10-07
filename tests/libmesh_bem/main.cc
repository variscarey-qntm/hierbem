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
 * @file main.cc
 * @brief Catch2 main function which initializes libMesh.
 */

#include <catch2/catch_session.hpp>
#include <libmesh/libmesh.h>

#include "test_common.h"

namespace HierBEM
{
  namespace LibMeshBEM
  {
    namespace testing
    {
      libMesh::LibMeshInit *libmesh_init = nullptr;
    }
  } // namespace LibMeshBEM
} // namespace HierBEM

int
main(int argc, char **argv)
{
  // libMesh only needs to see the program name; the remaining command line
  // arguments are interpreted by Catch2.
  int                  libmesh_argc = 1;
  libMesh::LibMeshInit init(libmesh_argc, argv);
  HierBEM::LibMeshBEM::testing::libmesh_init = &init;

  return Catch::Session().run(argc, argv);
}
