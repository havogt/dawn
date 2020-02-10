//---- Preprocessor defines ----
#define DAWN_GENERATED 1
#define DAWN_BACKEND_T CXXNAIVEICO
#include <driver-includes/unstructured_interface.hpp>

//---- Includes ----
#include "driver-includes/gridtools_includes.hpp"
using namespace gridtools::dawn;

//---- Globals ----

//---- Stencils ----
namespace dawn_generated {
namespace cxxnaiveico {
template <typename LibTag>
class fvm_nabla {
private:
  struct stencil_40 {
    dawn::mesh_t<LibTag> const& m_mesh;
    int m_k_size;
    dawn::edge_field_t<LibTag, double>& m_S_MXX;
    dawn::edge_field_t<LibTag, double>& m_zavgS_MXX;
    dawn::vertex_field_t<LibTag, double>& m_pp;
    dawn::vertex_field_t<LibTag, double>& m_pnabla_MXX;
    dawn::vertex_field_t<LibTag, double>& m_vol;
    dawn::sparse_vertex_field_t<LibTag, double>& m_sign;

  public:
    stencil_40(dawn::mesh_t<LibTag> const& mesh, int k_size,
               dawn::edge_field_t<LibTag, double>& S_MXX,
               dawn::edge_field_t<LibTag, double>& zavgS_MXX,
               dawn::vertex_field_t<LibTag, double>& pp,
               dawn::vertex_field_t<LibTag, double>& pnabla_MXX,
               dawn::vertex_field_t<LibTag, double>& vol,
               dawn::sparse_vertex_field_t<LibTag, double>& sign)
        : m_mesh(mesh), m_k_size(k_size), m_S_MXX(S_MXX), m_zavgS_MXX(zavgS_MXX), m_pp(pp),
          m_pnabla_MXX(pnabla_MXX), m_vol(vol), m_sign(sign) {}

    ~stencil_40() {}

    void sync_storages() {}

    void run() {
      using dawn::deref;
      ;
      {
        for(int k = 0 + 0; k <= (m_k_size == 0 ? 0 : (m_k_size - 1)) + 0 + 0; ++k) {
          for(auto const& loc : getCells(LibTag{}, m_mesh)) {
            int m_sparse_dimension_idx = 0;
            ::dawn::float_type __local_zavg_74 =
                ((::dawn::float_type)0.5 *
                 reduceVertexToEdge(LibTag{}, m_mesh, loc, (::dawn::float_type)0.0,
                                    [&](auto& lhs, auto const& red_loc) {
                                      lhs += m_pp(deref(LibTag{}, red_loc), k + 0);
                                      m_sparse_dimension_idx++;
                                      return lhs;
                                    }));
            m_zavgS_MXX(deref(LibTag{}, loc), k + 0) =
                (m_S_MXX(deref(LibTag{}, loc), k + 0) * __local_zavg_74);
            m_sparse_dimension_idx = 0;
            m_pnabla_MXX(deref(LibTag{}, loc), k + 0) = reduceEdgeToVertex(
                LibTag{}, m_mesh, loc, (::dawn::float_type)0.0,
                [&](auto& lhs, auto const& red_loc) {
                  lhs += (m_zavgS_MXX(deref(LibTag{}, red_loc), k + 0) *
                          m_sign(deref(LibTag{}, loc), m_sparse_dimension_idx, k + 0));
                  m_sparse_dimension_idx++;
                  return lhs;
                });
            m_pnabla_MXX(deref(LibTag{}, loc), k + 0) =
                (m_pnabla_MXX(deref(LibTag{}, loc), k + 0) / m_vol(deref(LibTag{}, loc), k + 0));
          }
        }
      }
      sync_storages();
    }
  };
  static constexpr const char* s_name = "fvm_nabla";
  stencil_40 m_stencil_40;

public:
  fvm_nabla(const fvm_nabla&) = delete;

  // Members

  fvm_nabla(const dawn::mesh_t<LibTag>& mesh, int k_size, dawn::edge_field_t<LibTag, double>& S_MXX,
            dawn::edge_field_t<LibTag, double>& S_MYY,
            dawn::edge_field_t<LibTag, double>& zavgS_MXX,
            dawn::edge_field_t<LibTag, double>& zavgS_MYY, dawn::vertex_field_t<LibTag, double>& pp,
            dawn::vertex_field_t<LibTag, double>& pnabla_MXX,
            dawn::vertex_field_t<LibTag, double>& pnabla_MYY,
            dawn::vertex_field_t<LibTag, double>& vol,
            dawn::sparse_vertex_field_t<LibTag, double>& sign)
      : m_stencil_40(mesh, k_size, S_MXX, zavgS_MXX, pp, pnabla_MXX, vol, sign) {}

  void run() {
    m_stencil_40.run();
    ;
  }
};
} // namespace cxxnaiveico
} // namespace dawn_generated
