#include <fstream>

#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/meshgenerator.h"
#include "atlas/option/Options.h"
#include "atlas/output/Gmsh.h"
#include "atlas_interface.hpp"

#include "fvm_nabla_hand_fixed.cpp"

int main() {
  // see atlas/numerics/fmv/Nabla.cc gradient_of_scalar

  // vol = dual_volumes
  // S = dual_normals
  // sign = see datastruct_module.f90 line 393ff

  // lonlat for initialization of the test field

  // octahedral: e.g. "O32"
  atlas::StructuredGrid structuredGrid = atlas::Grid("O32");
  atlas::StructuredMeshGenerator generator;
  auto mesh = generator.generate(structuredGrid);
  atlas::mesh::actions::build_edges(mesh);
  atlas::mesh::actions::build_node_to_edge_connectivity(mesh);
  atlas::mesh::actions::build_median_dual_mesh(mesh);

  int nb_levels = 1;
  atlas::functionspace::EdgeColumns fs_edges(mesh, atlas::option::levels(nb_levels));
  atlas::Field m_S_MXX{fs_edges.createField<double>(atlas::option::name("m_S_MXX"))};
  atlas::Field m_S_MYY{fs_edges.createField<double>(atlas::option::name("m_S_MYY"))};
  atlas::Field m_zavgS_MXX{fs_edges.createField<double>(atlas::option::name("m_zavgS_MXX"))};
  atlas::Field m_zavgS_MYY{fs_edges.createField<double>(atlas::option::name("m_zavgS_MYY"))};

  atlas::functionspace::NodeColumns fs_nodes(mesh, atlas::option::levels(nb_levels));
  atlas::Field m_pp{fs_nodes.createField<double>(atlas::option::name("m_pp"))};
  atlas::Field m_pnabla_MXX{fs_nodes.createField<double>(atlas::option::name("m_pnabla_MXX"))};
  atlas::Field m_pnabla_MYY{fs_nodes.createField<double>(atlas::option::name("m_pnabla_MYY"))};
  atlas::Field m_vol{fs_nodes.createField<double>(atlas::option::name("m_vol"))};
  constexpr int edges_per_node = 7; // TODO change for octahedral
  atlas::Field m_sign{fs_nodes.createField<double>(atlas::option::name("m_sign") |
                                                   atlas::option::variables(edges_per_node))};

  const auto vol_atlas = atlas::array::make_view<double, 1>(mesh.nodes().field("dual_volumes"));
  const auto S = atlas::array::make_view<double, 2>(mesh.edges().field("dual_normals"));

  constexpr int MXX = 0;
  constexpr int MYY = 1;

  {
    // all fields supported by dawn are 2 dimensional: (unstructured, lev)
    // S has dimensions (unstructured, [MMX/MMY])
    auto S_MXX = atlas::array::make_view<double, 2>(m_S_MXX);
    auto S_MYY = atlas::array::make_view<double, 2>(m_S_MYY);
    // TODO for nblevels > 1 we are missing a copy loop here
    for(int i = 0, size = mesh.edges().size(); i < size; ++i) {
      S_MXX(i, 0) = S(i, MXX);
      S_MYY(i, 0) = S(i, MYY);
    }
    auto vol = atlas::array::make_view<double, 2>(m_vol);
    // TODO for nblevels > 1 we are missing a copy loop here
    for(int i = 0, size = mesh.nodes().size(); i < size; ++i) {
      vol(i, 0) = vol_atlas(i);
    }

    // compute sign field
    auto node2edge_sign = atlas::array::make_view<double, 3>(m_sign);

    auto edge_flags = atlas::array::make_view<int, 1>(mesh.edges().flags());
    using Topology = atlas::mesh::Nodes::Topology;
    auto is_pole_edge = [&](size_t e) { return Topology::check(edge_flags(e), Topology::POLE); };

    for(std::size_t jnode = 0; jnode < mesh.nodes().size(); ++jnode) {
      auto const& node_edge_connectivity = mesh.nodes().edge_connectivity();
      auto const& edge_node_connectivity = mesh.edges().node_connectivity();
      for(std::size_t jedge = 0; jedge < node_edge_connectivity.cols(jnode); ++jedge) {
        auto iedge = node_edge_connectivity(jnode, jedge);
        auto ip1 = edge_node_connectivity(iedge, 0);
        if(jnode == ip1) {
          node2edge_sign(jnode, 0, jedge) = 1.;
        } else {
          node2edge_sign(jnode, 0, jedge) = -1.;
          if(is_pole_edge(iedge)) {
            node2edge_sign(jnode, 0, jedge) = 1.;
          }
        }
      }
    }

    auto pp = atlas::array::make_view<double, 2>(m_pp);
    // TODO for nblevels > 1 we are missing a copy loop here
    for(int i = 0, size = mesh.nodes().size(); i < size; ++i) {
      pp(i, 0) = 1;
    }
  }

  atlasInterface::Field<double> vol = atlas::array::make_view<double, 2>(m_vol);
  atlasInterface::Field<double> S_MXX = atlas::array::make_view<double, 2>(m_S_MXX);
  atlasInterface::Field<double> S_MYY = atlas::array::make_view<double, 2>(m_S_MYY);
  atlasInterface::Field<double> zavgS_MXX = atlas::array::make_view<double, 2>(m_zavgS_MXX);
  atlasInterface::Field<double> zavgS_MYY = atlas::array::make_view<double, 2>(m_zavgS_MYY);
  atlasInterface::Field<double> pp = atlas::array::make_view<double, 2>(m_pp);
  atlasInterface::Field<double> pnabla_MXX = atlas::array::make_view<double, 2>(m_pnabla_MXX);
  atlasInterface::Field<double> pnabla_MYY = atlas::array::make_view<double, 2>(m_pnabla_MYY);
  atlasInterface::SparseDimension<double> sign = atlas::array::make_view<double, 3>(m_sign);

  // TODO fill sign field

  atlas::output::Gmsh gmesh("mymesh.msh");
  gmesh.write(mesh);
  gmesh.write(m_pp);
  for(int i = 0; i < 1; ++i) {
    // in.metadata().set("step", i);

    dawn_generated::cxxnaiveico::fvm_nabla<atlasInterface::atlasTag>(
        mesh, nb_levels, S_MXX, S_MYY, zavgS_MXX, zavgS_MYY, pp, pnabla_MXX, pnabla_MYY, vol, sign)
        .run();
    gmesh.write(m_pnabla_MXX);
  }
}
