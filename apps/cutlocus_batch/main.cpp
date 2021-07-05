#include <map>
#include <stdio.h>
#include <thread>
#include <vector>
#include <yocto/yocto_cli.h>
#include <yocto/yocto_commonio.h>
#include <yocto/yocto_geometry.h>
using namespace std;
#include "libigl/include/igl/heat_geodesics.h"
#include <src/clio.h>
#include <src/cut_locus.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_opengl.h>
//#include <yocto_gui/yocto_window.h>
using namespace yocto;
void export_branches(const shape_data &mesh, const vector<node> &cut_graph) {
  if (cut_graph.size() > 0) {
    std::ofstream outfile;
    auto [branches, br_pos] = compute_branches(mesh, cut_graph);
    outfile.open("CL_On_Verts");
    for (auto i = 0; i < branches.size(); ++i) {
      outfile << cut_graph[branches[i].x].vid << " ";
      outfile << cut_graph[branches[i].y].vid << "\n";
    }
    outfile.close();
  }
}
void export_smoothed_branches(const shape_data &mesh,
                              const shape_topology &topology,
                              const vector<node> &cut_graph,
                              const vector<mesh_point> &new_ones,
                              const bool &draw) {
  auto node_pos = vector<vector<vec3f>>{};
  for (auto i = 0; i < cut_graph.size(); ++i) {
    if (cut_graph[i].is_cut || !cut_graph[i].is_in_restricted)
      continue;
    auto curr_n = tid_normal(mesh.triangles, mesh.positions, new_ones[i].face);
    for (auto h = 1; h < cut_graph[i].neighbors.size(); ++h) {

      auto child_entry = cut_graph[i].neighbors[h];
      if (cut_graph[child_entry].is_cut ||
          !cut_graph[child_entry].is_in_restricted)
        continue;
      node_pos.push_back({eval_position(mesh, new_ones[i]),
                          eval_position(mesh, new_ones[child_entry])});
    }
  }
  std::ofstream outfile0;
  outfile0.open("CL_Smoothed");
  for (auto i = 0; i < node_pos.size(); ++i) {
    for (auto j = 0; j < node_pos[i].size() - 1; ++j) {
      outfile0 << node_pos[i][j].x << " ";
      outfile0 << node_pos[i][j].y << " ";
      outfile0 << node_pos[i][j].z << " ";
      outfile0 << node_pos[i][j + 1].x << " ";
      outfile0 << node_pos[i][j + 1].y << " ";
      outfile0 << node_pos[i][j + 1].z << " \n";
    }
  }

  outfile0.close();
}
vector<float> compute_distance_field(const shape_data &mesh,
                                     const geodesic_solver &solver,
                                     const int source, const int method) {
  auto distances = vector<float>{};
  switch (method) {
  case VTP:
    distances = exact_geodesic_distance(mesh.triangles, mesh.positions, source);
    break;
  case heat_method: {
    auto heat_solver = igl::HeatGeodesicsData<double>{};
    auto [V, F] = libigl_wrapper(mesh.positions, mesh.triangles);
    igl::heat_geodesics_precompute(V, F, heat_solver);
    distances = heat_geodesic(heat_solver, source);

  } break;
  case graph_solver:
    distances = compute_geodesic_distances(solver, {source});
    break;
  }
  return distances;
}
void compute_cut_locus(const shape_data &mesh, const shape_topology &topology,
                       const geodesic_solver &solver,
                       const dual_geodesic_solver &dual_solver,
                       const discr_diff_op &operators, const int source,
                       const int method, const bool use_grad,
                       const bool use_lap, const float &angle_threshold,
                       const float &lap_threshold, const int depth_threshold,
                       const bool growing, const bool smoothing,
                       const int smoothing_grad_it,
                       const int smoothing_normal_it, vector<node> &cut_graph,
                       vector<mesh_point> &smoothed) {
  auto scaling_factor = 1.f;
  auto distances = compute_distance_field(mesh, solver, source, method);
  auto gradients = compute_grad(solver, mesh.triangles, mesh.positions,
                                mesh.normals, operators.Grad, distances);
  update_angle_deviation(mesh, topology, solver, gradients, cut_graph);
  auto [lap, entry_of_min, max_value] =
      compute_laplacian_with_extrema(operators, distances);
  auto max_of_dist = field_maximum(distances);
  lap[max_of_dist] = {flt_min, max_of_dist};
  spanning_tree(mesh, topology, solver, lap, gradients, max_of_dist, cut_graph);
  update_restricted(angle_threshold, lap_threshold, use_grad, use_lap, growing,
                    max_of_dist, cut_graph);

  auto genus =
      (int)(-(mesh.positions.size() - mesh.triangles.size() / 2) / 2 + 1);
  if (genus > 0) {
    auto source_point =
        make_point_from_vert(mesh.triangles, topology.v2t, source);
    auto edge_mask = vector<vec3i>(mesh.triangles.size());
    auto mixed_tree = vector<node>{};
    auto cotree = vector<node>{};
    mixed_spanning_tree(mesh, topology, solver, cut_graph, distances, lap,
                        max_value, max_of_dist, mixed_tree, edge_mask);
    spanning_cotree(mesh, topology, distances, source_point, cotree, edge_mask);
    auto generators = vector<vec2i>{};
    vector<vector<int>> links = {};
    for (auto i = 0; i < edge_mask.size(); ++i) {
      for (auto j = 0; j < 3; ++j) {
        if (edge_mask[i][j] == 0) {
          generators.push_back(
              vec2i{mesh.triangles[i][j], mesh.triangles[i][(j + 1) % 3]});
          auto nei = topology.adjacencies[i][j];
          auto k = find(topology.adjacencies[nei], i);
          edge_mask[nei][k] = 2;
        }
      }
    }
    links = find_generators(mesh, topology, solver, mixed_tree, cotree,
                            generators, true);
    update_tree(links, lap, cut_graph);

    update_restricted(angle_threshold, lap_threshold, use_grad, use_lap,
                      growing, max_of_dist, cut_graph);
  }
  if (depth_threshold > 0)
    cut_short_branches_in_restricted(cut_graph, depth_threshold);

  if (smoothing) {
    smoothed =
        smoothed_cut_locus(mesh, topology, solver, dual_solver, cut_graph,
                           gradients, smoothing_grad_it, scaling_factor, true);
    smoothed_cut_locus(mesh, topology, solver, dual_solver, cut_graph,
                       gradients, smoothing_normal_it, scaling_factor, false,
                       smoothed);
  }
}
int get_method(const std::string s) {
  if (s == "VTP")
    return VTP;
  else if (s == "Heat")
    return heat_method;
  else if (s == "Graph")
    return graph_solver;
  else
    return -1;
}
int main(int num_args, const char *args[]) {
  string filename = "data/mesh.ply", source = "0", method = "VTP";
  auto use_grad = false, use_lap = false, growing = false, cut = false,
       smoothing = false;
  auto angle_threshold = 0.f, lap_threshold = 0.f;
  auto depth_threshold = 0, smoothing_grad_it = 0, smoothing_normal_it = 0;
  auto cli = make_cli("Cut Locus", "Batch Program for Computing the cut locus");
  add_argument(cli, "mesh", filename, "Model name");
  add_argument(cli, "source", source, "Vert index");
  add_argument(cli, "method", method, "Mehtod for geodesic distances");
  add_option(cli, "use_grad", use_grad, "Use gradient filtering");
  add_option(cli, "use_lap", use_lap, "Use Laplacian filtering");
  add_option(cli, "cut", cut, "Cut short branches");
  add_option(cli, "growing", growing, "Set growing policy");
  add_option(cli, "smoothed", smoothing, "Smoothed CL");
  parse_cli_and_handle_errors(cli, vector<string>{args, args + num_args});
  auto mesh = shape_data{};
  auto topology = shape_topology{};
  auto operators = discr_diff_op{};
  auto solver = geodesic_solver{};
  auto dual_solver = dual_geodesic_solver{};
  auto error = string{};

  if (!load_mesh(filename, mesh, topology, operators, solver, dual_solver,
                 error))
    exit(1);

  if (use_grad) {
    std::cout << "\n"
              << "Choose the threshold for the gradient filtering: ";
    std::cin >> angle_threshold;
  }
  if (use_lap) {
    std::cout << "\n"
              << "Choose the threshold for the laplacian filtering: ";
    std::cin >> lap_threshold;
  }
  if (cut) {
    std::cout << "\n"
              << "Choose the minimum length for branches: ";
    std::cin >> depth_threshold;
  }
  if (smoothing) {
    std::cout << "\n"
              << "Choose the number of smoothing iterations using gradient: ";
    std::cin >> smoothing_grad_it;
    std::cout << "\n"
              << "Choose the number of smoothing iterations using normal: ";
    std::cin >> smoothing_normal_it;
  }
  auto cut_graph = vector<node>{};
  auto smoothed = vector<mesh_point>{};
  compute_cut_locus(mesh, topology, solver, dual_solver, operators,
                    stoi(source), get_method(method), use_grad, use_lap,
                    angle_threshold, lap_threshold, depth_threshold, growing,
                    smoothing, smoothing_grad_it, smoothing_normal_it,
                    cut_graph, smoothed);
  export_branches(mesh, cut_graph);
  if (smoothing)
    export_smoothed_branches(mesh, topology, cut_graph, smoothed, false);
}