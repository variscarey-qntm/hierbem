# Copyright (C) 2026 HierBEM contributors
#
# This file is part of the HierBEM library.
#
# HierBEM is free software: you can use it, redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version. The full text of the license can be found in the
# file LICENSE at the top level directory of HierBEM.

# This example is only built with the libMesh backend (-DHBEM_USE_LIBMESH=ON).

set(_target solve-electrostatic-plate)

add_executable(${_target} solve-electrostatic-plate.cc)

target_compile_definitions(${_target}
                           PRIVATE SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(${_target} hierbem_libmesh)
set_target_properties(${_target} PROPERTIES BUILD_RPATH
                                            "${LIBMESH_LIBRARY_DIRS}")
