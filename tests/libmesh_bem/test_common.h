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
 * @file test_common.h
 * @brief Shared helpers for the libMesh backend tests.
 */

#ifndef HIERBEM_TESTS_LIBMESH_BEM_TEST_COMMON_H_
#define HIERBEM_TESTS_LIBMESH_BEM_TEST_COMMON_H_

#include <libmesh/libmesh.h>

namespace HierBEM
{
  namespace LibMeshBEM
  {
    namespace testing
    {
      extern libMesh::LibMeshInit *libmesh_init;
    }
  } // namespace LibMeshBEM
} // namespace HierBEM

#endif // HIERBEM_TESTS_LIBMESH_BEM_TEST_COMMON_H_
