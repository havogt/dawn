#include "dawn/CodeGen/CXXNaive-ico/CXXNaiveCodeGen.h"
#include "dawn/CodeGen/CodeGen.h"
#include <fstream>
#include <sstream>
#include <string>

#include "code_generator.hpp"

template <typename CG>
void dump(std::string filename, dawn::codegen::stencilInstantiationContext const& ctx) {
  std::ofstream of(filename);
  dawn::DiagnosticsEngine diagnostics;
  CG generator(ctx, diagnostics, 0);
  auto tu = generator.generateCode();

  std::ostringstream ss;
  for(auto const& macroDefine : tu->getPPDefines())
    ss << macroDefine << "\n";

  ss << tu->getGlobals();
  for(auto const& s : tu->getStencils())
    ss << s.second;
  of << ss.str();
  of.close();
}

int main() {
  for(auto&& stencil_instantiation_context :
      dawn::unstructured_codegen_driver::getCodegenRegistry()) {
    auto&& [name, si] = *(stencil_instantiation_context.begin());
    dump<dawn::codegen::cxxnaiveico::CXXNaiveIcoCodeGen>("generated_" + name + ".hpp",
                                                         stencil_instantiation_context);
  }
}
