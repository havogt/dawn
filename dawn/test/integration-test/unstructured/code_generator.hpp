#pragma once

#include "dawn/CodeGen/CodeGen.h"
#include "dawn/IIR/Stencil.h"
#include <string>
#include <vector>

#define DAWN_UNSTRUCTURED_CODEGEN(DAWN_UNSTRUCTURED_CODEGEN_NAME)                                  \
  class DAWN_UNSTRUCTURED_CODEGEN_NAME {                                                           \
  private:                                                                                         \
    static dawn::unstructured_codegen_driver::internal::CodeGenInfo const code_gen_info_           \
        __attribute__((unused));                                                                   \
    dawn::codegen::stencilInstantiationContext generate();                                         \
  };                                                                                               \
  dawn::unstructured_codegen_driver::internal::CodeGenInfo const                                   \
      DAWN_UNSTRUCTURED_CODEGEN_NAME::code_gen_info_ =                                             \
          dawn::unstructured_codegen_driver::internal::registerCodegen(                            \
              DAWN_UNSTRUCTURED_CODEGEN_NAME{}.generate());                                        \
  dawn::codegen::stencilInstantiationContext DAWN_UNSTRUCTURED_CODEGEN_NAME::generate()

namespace dawn::unstructured_codegen_driver {

// Global map of name to stencilInstantiationContext.
// TODO move definition to .cpp
std::vector<dawn::codegen::stencilInstantiationContext>& getCodegenRegistry();

namespace internal {
struct CodeGenInfo {};

// Register a stencilInstantiationContext in the global map.
CodeGenInfo registerCodegen(dawn::codegen::stencilInstantiationContext&& si_context);
} // namespace internal

} // namespace dawn::unstructured_codegen_driver
