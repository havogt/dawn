#include "dawn/CodeGen/CXXNaive-ico/CXXNaiveCodeGen.h"
#include "dawn/CodeGen/CXXNaive/CXXNaiveCodeGen.h"
#include "dawn/Unittest/IIRBuilder.h"
#include <string>

#include "code_generator.hpp"
#include "gtest/gtest.h"

namespace {

DAWN_UNSTRUCTURED_CODEGEN(CopyCell) {
  using namespace dawn::iir;
  using LocType = dawn::ast::LocationType;

  UnstructuredIIRBuilder b;
  auto in_f = b.field("in_field", LocType::Cells);
  auto out_f = b.field("out_field", LocType::Cells);

  return b.build(b.stencil(
      b.multistage(LoopOrderKind::Parallel,
                   b.stage(b.doMethod(dawn::sir::Interval::Start, dawn::sir::Interval::End,
                                      b.stmt(b.assignExpr(b.at(out_f), b.at(in_f))))))));
}

// TEST(AtlasIntegrationTestCompareOutput, CopyCell) {
//   // Setup a 32 by 32 grid of quads and generate a mesh out of it
//   atlas::StructuredGrid structuredGrid = atlas::Grid("L32x32");
//   atlas::StructuredMeshGenerator generator;
//   auto mesh = generator.generate(structuredGrid);

//   // We only need one vertical level
//   size_t nb_levels = 1;

//   // Create input (on cells) and output (on cells) fields
//   atlas::functionspace::CellColumns fs_cells(mesh, atlas::option::levels(nb_levels));
//   atlas::Field in{fs_cells.createField<double>(atlas::option::name("in"))};
//   atlas::Field out{fs_cells.createField<double>(atlas::option::name("out"))};

//   // Make views on the fields (needed to access the field like an array)
//   atlasInterface::Field<double> in_v = atlas::array::make_view<double, 2>(in);
//   atlasInterface::Field<double> out_v = atlas::array::make_view<double, 2>(out);

//   // Initialize fields with data
//   for(int cell_idx = 0; cell_idx < mesh.cells().size(); ++cell_idx) {
//     in_v(cell_idx, 0) = 1.0;
//     out_v(cell_idx, 0) = -1.0;
//   }

//   // Run the stencil
//   dawn_generated::cxxnaiveico::copyCell<atlasInterface::atlasTag>(mesh,
//   static_cast<int>(nb_levels),
//                                                                   in_v, out_v)
//       .run();

//   // Check correctness of the output
//   for(int cell_idx = 0; cell_idx < mesh.cells().size(); ++cell_idx)
//     ASSERT_EQ(out_v(cell_idx, 0), 1.0);
// }
} // namespace
