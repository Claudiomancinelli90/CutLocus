#pragma once
#include "geometry.h"
#include "libigl/include/igl/heat_geodesics.h"
// #include <Eigen/Dense>
// #include <Eigen/Sparse>
#include <deque>
#include <iostream>
#include <src/logging.h>
#include <yocto/yocto_color.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_shape.h>

// TODO(giacomo): is this necessary?
inline bool operator==(const mesh_point &a, const mesh_point &b) {
  return (a.face == b.face) && (a.uv == b.uv);
}

inline bool check_point(const mesh_point &point) {
  assert(point.face != -1);
  assert(point.uv.x >= 0);
  assert(point.uv.y >= 0);
  assert(point.uv.x <= 1);
  assert(point.uv.y <= 1);
  return true;
}
struct node {
  int vid = -1;    // id of the vertex
  float lap = 0.f; // laplacian at vertex
  float angle_deviation = -100;
  vector<int> neighbors = {}; // children
  bool is_in_restricted = false;
  bool is_cut = false;
  bool unprunable = false;
};

// struct bezier_mesh {
//   vector<vec3i> triangles = {};
//   vector<vec3i> adjacencies = {};
//   vector<vec3f> positions = {};
//   vector<vec3f> normals = {};
//   vector<vec2f> texcoords = {};
//   // Additional data for experimental algorithms
//   geodesic_solver solver = {};
//   vector<vector<int>> v2t = {};
//   vector<vector<float>> angles = {};
//   vector<float> total_angles = {};
//   dual_geodesic_solver dual_solver = {};
//   Eigen::SparseMatrix<double, 1> Grad;
//   Eigen::SparseMatrix<double> Lap;
//   vector<Eigen::MatrixXd> a;
//   vector<vector<vector<pair<int, float>>>> ring;
//   // vector<vector<vec2f>> christoffel;
//   vector<vector<vector<float>>> christoffel;
//   float avg_edge_length = 0.f;
// };

enum type_of_scalar_field { dist, lap };
enum type_of_distance_field { VTP, heat_method, graph_solver };

vector<float> heat_geodesic(const igl::HeatGeodesicsData<double> &solver,
                            const int &source);

int field_maximum(const vector<float> &distances);

geodesic_path compute_geodesic_path(const shape_data &mesh,
                                    const shape_topology &topology,
                                    const dual_geodesic_solver &dual_solver,
                                    const geodesic_solver &solver,
                                    const mesh_point &start,
                                    const mesh_point &end);
std::tuple<vector<pair<float, int>>, int, float>
compute_laplacian_with_extrema(const discr_diff_op &operators,
                               const vector<float> &distances);

void update_angle_deviation(const shape_data &mesh,
                            const shape_topology &topology,
                            const geodesic_solver &solver,
                            const vector<vec3f> &gradients,
                            vector<node> &graph);

void spanning_tree(const shape_data &mesh, const shape_topology &topology,
                   const geodesic_solver &solver,
                   const vector<pair<float, int>> &lap,
                   const vector<vec3f> &gradients, const int entry_of_root,
                   vector<node> &graph);
void mixed_spanning_tree(const shape_data &mesh, const shape_topology &topology,
                         const geodesic_solver &solver,
                         const vector<node> &lap_graph,
                         const vector<float> &distances,
                         const vector<pair<float, int>> &lap,
                         const float &max_lap, const int entry_of_min,
                         vector<node> &mixed_graph, vector<vec3i> &edge_mask);

void spanning_cotree(const shape_data &mesh, const shape_topology &topology,
                     const vector<float> &distances, const mesh_point &source,
                     vector<node> &graph, vector<vec3i> &edge_mask);
vector<vector<int>>
find_generators(const shape_data &mesh, const shape_topology &topology,
                const geodesic_solver &solver,
                const vector<node> &cut_graph_mixed, const vector<node> &cotree,
                const vector<vec2i> &generators, const bool fast);

void update_tree(const vector<vector<int>> &links,
                 const vector<pair<float, int>> &lap, vector<node> &cut_tree);

void cut_short_branches(vector<node> &curr_tree, const int depth);
void cut_short_branches_in_restricted(vector<node> &curr_tree, const int depth);

void update_restricted(const float &grad_threshold, const float &lap_threshold,
                       const bool &use_grad, const bool &use_lap,
                       const bool &propagate, const int entry_of_root,
                       vector<node> &cut_tree);

void uncut_tree(vector<node> &curr_tree);
void uncut_tree_in_restricted(vector<node> &curr_tree);
void thinning(const shape_data &mesh, const geodesic_solver &solver,
              const vector<pair<float, int>> &lap, vector<node> &cut_graph);
void thinning_in_restricted(const shape_data &mesh,
                            const geodesic_solver &solver,
                            const vector<pair<float, int>> &lap,
                            vector<node> &cut_graph);

vector<mesh_point>
smoothed_cut_locus(const shape_data &mesh, const shape_topology &topology,
                   const geodesic_solver &solver,
                   const dual_geodesic_solver &dual_solver,
                   vector<node> &cut_graph, const vector<vec3f> &gradients,
                   const int it, float &factor, const bool &use_gradient);
void smoothed_cut_locus(const shape_data &mesh, const shape_topology &topology,
                        const geodesic_solver &solver,
                        const dual_geodesic_solver &dual_solver,
                        vector<node> &cut_graph, const vector<vec3f> &gradients,
                        const int it, float &factor, const bool &use_gradient,
                        vector<mesh_point> &mapping);

std::tuple<vector<vec2i>, vector<vec3f>, vector<vec2i>, vector<vec3f>>
compute_splitted_branches(const shape_data &mesh,
                          const vector<node> &cut_graph);
std::pair<vector<vec2i>, vector<vec3f>>
compute_branches(const shape_data &mesh, const vector<node> &cut_graph);
std::pair<vector<vec2i>, vector<vec3f>>
compute_cobranches(const shape_data &mesh, const shape_topology &topology,
                   const vector<node> &cut_graph);

std::pair<vector<vec2i>, vector<vec3f>>
compute_smoothed_branches(const shape_data &mesh,
                          const shape_topology &topology,
                          const vector<node> &cut_graph,
                          const vector<mesh_point> &new_one, const bool &draw);

vec3f polar_basis(const geodesic_solver &solver, const vector<vec3f> &positions,
                  const vector<vec3f> &normals, int vid);

vector<vec2i> export_cobranches(const shape_data &mesh,
                                const vector<node> &cut_graph);
inline vec3f make_bary(const vec2f &uv) {
  return vec3f{1 - uv.x - uv.y, uv.x, uv.y};
}

inline vec3f eval_position(const shape_data &mesh, const mesh_point &point) {
  return eval_position(mesh.triangles, mesh.positions, point);
}

inline vector<vec3f> path_positions(const shape_data &mesh,
                                    const shape_topology &topology,
                                    const geodesic_path &path) {
  return path_positions(path, mesh.triangles, mesh.positions,
                        topology.adjacencies);
}

inline void parallel_transport(const shape_data &mesh,
                               const shape_topology &topology,
                               const geodesic_solver &solver, vec3f &v,
                               const int from, const int to, const int mode) {
  parallel_transp(solver, topology.angles, topology.total_angles,
                  mesh.triangles, mesh.positions, topology.adjacencies, v,
                  mesh.normals, from, to, mode);
}
