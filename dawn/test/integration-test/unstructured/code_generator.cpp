#include "code_generator.hpp"

namespace dawn::unstructured_codegen_driver {

// Global map of name to stencilInstantiationContext.
// TODO move definition to .cpp
std::vector<dawn::codegen::stencilInstantiationContext>& getCodegenRegistry() {
  static std::vector<dawn::codegen::stencilInstantiationContext> ctxs;
  return ctxs;
}

namespace internal {
// Register a stencilInstantiationContext in the global map.
CodeGenInfo registerCodegen(dawn::codegen::stencilInstantiationContext&& ctx) {
  getCodegenRegistry().emplace_back(ctx);
  return {};
}
} // namespace internal

} // namespace dawn::unstructured_codegen_driver
