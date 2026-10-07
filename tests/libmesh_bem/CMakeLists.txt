# Copyright (C) 2026 HierBEM contributors
#
# This file is part of the HierBEM library.
#
# HierBEM is free software: you can use it, redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version. The full text of the license can be found in the
# file LICENSE at the top level directory of HierBEM.

# Unit tests for the libMesh backend (-DHBEM_USE_LIBMESH=ON).

include(FetchContent)
FetchContent_Declare(Catch2 SOURCE_DIR
                            "${CMAKE_SOURCE_DIR}/extern/Catch2-3.4.0/")
FetchContent_MakeAvailable(Catch2)

set(_target libmesh-bem-tests)

add_executable(
  ${_target} main.cc test_sauter_quadrature_triangle.cc test_gmsh_io.cc
             test_electrostatic_plate.cc test_lagrange_assembly.cc
             test_single_layer_entry.cc)

target_compile_definitions(
  ${_target}
  PRIVATE TEST_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
          EXAMPLE_MESH="${CMAKE_SOURCE_DIR}/examples/libmesh-electrostatic-plate/unit_square.msh"
)

target_link_libraries(${_target} hierbem_libmesh Catch2::Catch2)
set_target_properties(${_target} PROPERTIES BUILD_RPATH
                                            "${LIBMESH_LIBRARY_DIRS}")

add_test(NAME ${_target} COMMAND ${_target})
