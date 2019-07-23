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

yoda_add_library(
  NAME DawnSerializer
  SOURCES ASTSerializer.h
          ASTSerializer.cpp
          IIRSerializer.h
          IIRSerializer.cpp
          SIRSerializer.h
          SIRSerializer.cpp
  OBJECT
)

yoda_combine_libraries(
  NAME DawnSerializer
  OBJECTS DawnSerializerObjects
  DEPENDS DawnIIRStatic DawnOptimizerStatic
  INSTALL_DESTINATION ${DAWN_INSTALL_LIB_DIR}
  VERSION ${DAWN_VERSION}
  EXPORT_GROUP DawnTargets
)

add_dependencies(DawnSerializerObjects IIR_Proto_Generated)
target_include_directories(DawnSerializerObjects PUBLIC $<TARGET_PROPERTY:DawnIIRStatic,INCLUDE_DIRECTORIES>)

target_include_directories(DawnSerializerObjects PUBLIC /home/vogtha/projects/stencils_mlir/setup/llvm-project/install/include)

target_link_directories(DawnSerializerShared PUBLIC /home/vogtha/projects/stencils_mlir/setup/llvm-project/install/lib)
target_link_libraries(DawnSerializerShared PUBLIC MLIRIR MLIRStencil LLVMCore LLVMSupport MLIRSupport)
target_link_libraries(DawnSerializerShared PUBLIC z rt dl tinfo pthread m xml2)

target_link_directories(DawnSerializerStatic PUBLIC /home/vogtha/projects/stencils_mlir/setup/llvm-project/install/lib)
target_link_libraries(DawnSerializerStatic PUBLIC MLIRIR MLIRStencil LLVMCore LLVMSupport MLIRSupport)
target_link_libraries(DawnSerializerStatic PUBLIC z rt dl tinfo pthread m xml2)



#target_include_directories(DawnSerializerObjects PUBLIC ${sir_proto_include_dirs})
install(EXPORT DawnTargets NAMESPACE Dawn:: DESTINATION ${DAWN_INSTALL_CMAKE_DIR})
