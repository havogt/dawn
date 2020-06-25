//===--------------------------------------------------------------------------------*- C++ -*-===//
//                         _       _
//                        | |     | |
//                    __ _| |_ ___| | __ _ _ __   __ _
//                   / _` | __/ __| |/ _` | '_ \ / _` |
//                  | (_| | || (__| | (_| | | | | (_| |
//                   \__, |\__\___|_|\__,_|_| |_|\__, | - GridTools Clang DSL
//                    __/ |                       __/ |
//                   |___/                       |___/
//
//
//  This file is distributed under the MIT License (MIT).
//  See LICENSE.txt for details.
//
//===------------------------------------------------------------------------------------------===//
#include <gridtools/stencil_composition/sid/concept.hpp>
#define DAWN_GENERATED 1
#define GRIDTOOLS_DAWN_HALO_EXTENT 3
#define GT_VECTOR_LIMIT_SIZE 30

#undef FUSION_MAX_VECTOR_SIZE
#undef FUSION_MAX_MAP_SIZE
#define FUSION_MAX_VECTOR_SIZE GT_VECTOR_LIMIT_SIZE
#define FUSION_MAX_MAP_SIZE FUSION_MAX_VECTOR_SIZE
#define BOOST_MPL_LIMIT_VECTOR_SIZE FUSION_MAX_VECTOR_SIZE
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS

#include <gtest/gtest.h>
#include "test/integration-test/CodeGen/Macros.hpp"
#include "driver-includes/verify.hpp"
#include "test/integration-test/CodeGen/Options.hpp"
#include "test/integration-test/CodeGen/generated/copy_stencil_c++-naive.cpp"

#ifndef OPTBACKEND
#define OPTBACKEND gt
#endif

// clang-format off
#include INCLUDE_FILE(test/integration-test/CodeGen/generated/copy_stencil_,OPTBACKEND.cpp)
// clang-format on

using namespace dawn;
TEST(copy_stencil, test) {
  domain dom(Options::getInstance().m_size[0], Options::getInstance().m_size[1],
             Options::getInstance().m_size[2]);
  dom.set_halos(0, 0, 0, 0, 0, 0);
  // dom.set_halos(halo::value, halo::value, halo::value, halo::value, 0, 0);

  verifier verif(dom);

  meta_data_t meta_data(dom.isize(), dom.jsize(), dom.ksize() /* + 1 */);
  storage_t in(meta_data, "in"), out_gt(meta_data, "out-gt"), out_naive(meta_data, "out-naive");

  verif.fillMath(8.0, 2.0, 1.5, 1.5, 2.0, 4.0, in);
  verif.fill(-1.0, out_gt, out_naive);

  // dawn_generated::OPTBACKEND::copy_stencil copy_gt(dom);
  auto copy_gt = dawn_generated::OPTBACKEND::make_copy_stencil(dom);
  auto copy_naive = dawn_generated::cxxnaive::make_copy_stencil(dom);

  copy_gt(in, out_gt);
  copy_naive(in, out_naive);

  ASSERT_TRUE(verif.verify(in, out_naive));
  ASSERT_TRUE(verif.verify(out_gt, out_naive));
}

// TEST(copy_stencil, test_with_c_array) {
//   domain dom(Options::getInstance().m_size[0], Options::getInstance().m_size[1],
//              Options::getInstance().m_size[2]);
//   dom.set_halos(halo::value, halo::value, halo::value, halo::value, 0, 0);

//   // meta_data_t meta_data(dom.isize(), dom.jsize(), dom.ksize() + 1);
//   // storage_t in(meta_data, "in"), out_gt(meta_data, "out-gt"), out_naive(meta_data, "out-naive");

//   constexpr int size = 10;

//   float in[size][size][size];
//   float out_gt[size][size][size];
//   float out_naive[size][size][size];

//   for(std::size_t i = 0; i < size; ++i)
//     for(std::size_t j = 0; j < size; ++j)
//       for(std::size_t k = 0; k < size; ++k) {
//         in[i][j][k] = i + j * 10 + k * 100;
//         out_gt[i][j][k] = -1.;
//         out_naive[i][j][k] = -1.;
//       }

//   static_assert(gridtools::sid::is_sid<decltype(in)>{});

//   // // dawn_generated::OPTBACKEND::copy_stencil copy_gt(dom);
//   // auto copy_gt = dawn_generated::OPTBACKEND::make_copy_stencil<decltype(in), decltype(out_gt)>(dom);
//   auto copy_naive = dawn_generated::cxxnaive::make_copy_stencil(dom);

//   // copy_gt(in, out_gt);
//   copy_naive(in, out_naive);

//   // for(std::size_t i = 0; i < size; ++i)
//   //   for(std::size_t j = 0; j < size; ++j)
//   //     for(std::size_t k = 0; k < size; ++k) {
//   //       ASSERT_EQ(out_naive[i][j][k], out_gt[i][j][k]);
//   //     }
// }
