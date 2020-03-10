##===------------------------------------------------------------------------------*- CMake -*-===##
##                          _
##                         | |
##                       __| | __ ___      ___ ___
##                      / _` |/ _` \ \ /\ / / '_  |
##                     | (_| | (_| |\ V  V /| | | |
##                      \__,_|\__,_| \_/\_/ |_| |_| - Compiler Toolchain
##
##
##  This file is distributed under the MIT License (MIT).
##  See LICENSE.txt for details.
##
##===------------------------------------------------------------------------------------------===##

find_program(SH_PROGRAM sh)

function(generate_from_dawn4py)
  set(options)
  set(oneValueArgs NAME)
  set(multiValueArgs)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  add_custom_target(
    generate_${ARG_NAME}
    ALL
    COMMAND ${CMAKE_COMMAND} -E env
    PYTHONPATH=${DAWN4PY_MODULE_DIR}:${PROTOBUF_PYTHON_DIR}
    ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_NAME}.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}
  )
endfunction()

generate_from_dawn4py(NAME fvm_nabla)
add_executable(fvm_nabla_atlas_driver fvm_nabla_atlas_driver.cpp)
target_include_directories(fvm_nabla_atlas_driver PRIVATE ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY})
target_include_directories(fvm_nabla_atlas_driver PRIVATE ${CMAKE_SOURCE_DIR}/src)

target_add_dawn_standard_props(fvm_nabla_atlas_driver)
target_link_libraries(fvm_nabla_atlas_driver ${PROJECT_NAME} eckit atlas)
# target_link_libraries(fvm_nabla_atlas_drive GridTools::gridtools) # TODO shouldn't be needed
