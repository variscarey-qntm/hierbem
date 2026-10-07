# Copyright (C) 2026 HierBEM contributors
#
# This file is part of the HierBEM library.
#
# HierBEM is free software: you can use it, redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version. The full text of the license can be found in the
# file LICENSE at the top level directory of HierBEM.
# ------------------------------------------------------------------------------
#
# Build configuration for the libMesh backend of HierBEM, which provides
# Galerkin BEM on triangular surface meshes. It is enabled with
# -DHBEM_USE_LIBMESH=ON and does not depend on deal.II or CUDA.
#
# libMesh is located via pkg-config. Pass -DLIBMESH_DIR=<libmesh_prefix> or set
# the environment variable LIBMESH_DIR if libMesh is not installed in a
# standard location. The libMesh build method (opt, devel, dbg, ...) is
# selected with -DHBEM_LIBMESH_METHOD=<method> (default: opt).

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(HBEM_LIBMESH_METHOD
    "opt"
    CACHE STRING "libMesh build method (opt, devel, dbg, prof, oprof)")

if(NOT LIBMESH_DIR AND DEFINED ENV{LIBMESH_DIR})
  set(LIBMESH_DIR $ENV{LIBMESH_DIR})
endif()
if(LIBMESH_DIR)
  list(PREPEND CMAKE_PREFIX_PATH ${LIBMESH_DIR})
endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(LIBMESH REQUIRED IMPORTED_TARGET
                  libmesh-${HBEM_LIBMESH_METHOD})
message(STATUS "Found libMesh ${LIBMESH_VERSION} (${HBEM_LIBMESH_METHOD})")

# libMesh may be built with MPI. In that case, compile with the same MPI
# headers and libraries.
find_package(MPI COMPONENTS CXX)

add_library(
  hierbem_libmesh SHARED
  ${CMAKE_SOURCE_DIR}/src/libmesh_bem/sauter_quadrature_triangle.cc
  ${CMAKE_SOURCE_DIR}/src/libmesh_bem/triangle_pair.cc
  ${CMAKE_SOURCE_DIR}/src/libmesh_bem/galerkin_assembly.cc
  ${CMAKE_SOURCE_DIR}/src/libmesh_bem/gmsh_io.cc)
target_include_directories(hierbem_libmesh
                           PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(hierbem_libmesh PUBLIC PkgConfig::LIBMESH)
if(MPI_CXX_FOUND)
  target_link_libraries(hierbem_libmesh PUBLIC MPI::MPI_CXX)
endif()
target_compile_options(hierbem_libmesh PRIVATE -Wall -Wextra)
set_target_properties(hierbem_libmesh PROPERTIES BUILD_RPATH
                                                 "${LIBMESH_LIBRARY_DIRS}")

add_subdirectory(${CMAKE_SOURCE_DIR}/examples/libmesh-electrostatic-plate
                 ${CMAKE_BINARY_DIR}/examples/libmesh-electrostatic-plate)

enable_testing()
add_subdirectory(${CMAKE_SOURCE_DIR}/tests/libmesh_bem
                 ${CMAKE_BINARY_DIR}/tests/libmesh_bem)
