#pragma once
// #include <Eigen/Dense>
// #include <Eigen/Sparse>
#include <chrono>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::pair;
using std::string;
using std::vector;

} // namespace yocto
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <VTP/geodesic_algorithm_exact.h>
#include <VTP/geodesic_mesh.h>
#include <iostream>
#include <stdio.h>
#include <unordered_set>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_shape.h>
using namespace yocto;
enum transport_mode { V2V, V2T, T2V, T2T };
enum simplex { vert, edge };
enum field { exact, graph };

struct shape_topology {
  vector<vec3i> adjacencies = {};
  vector<vector<int>> v2t = {};
  vector<vector<float>> angles = {};
  vector<float> total_angles = {};
  float avg_edge_length = 0.f;
};

struct discr_diff_op {
  Eigen::SparseMatrix<double, 1> Grad;
  Eigen::SparseMatrix<double> Lap;
  vector<Eigen::MatrixXd> a;
  vector<vector<vector<pair<int, float>>>> ring;
  vector<vector<vec2f>> christoffel;
};
void clean_bary(mesh_point &sample);
vec3f tid_centroid(const vector<vec3i> &triangles,
                   const vector<vec3f> &positions, const int pid);

vec3f tri_bary_coords(const vec3f &v0, const vec3f &v1, const vec3f &v2,
                      const vec3f &p);
vec3f project_vec(const vec3f &v, const vec3f &n);

mesh_point make_point_from_vert(const vector<vec3i> &triangles,
                                const vector<vector<int>> &v2t, const int vid);
int forced_vert_from_point(const vector<vec3i> &triangles, const mesh_point &p);
vector<float> subdivide_angles(const int number_of_subdivision);

vec3f rot_vect(const vec3f &p, const vec3f &axis, const float angle);

mesh_point make_mesh_point(const vector<vec3i> &triangles,
                           const vector<vector<int>> &v2t, const int vid,
                           const int tid = -1);

float angle_in_tangent_space(const geodesic_solver &solver,
                             const vector<vec3f> &positions, const vec3f &v,
                             const int vid, const vec3f &n);
// vector<int> k_ring(const geodesic_solver &solver, const int vid, const int k,
//                    bool mesh_connectivity);
// vector<int> k_ring(const geodesic_solver &solver, const int vid,
//                    const int vid_on_boundary);
// vector<int> k_ring(const geodesic_solver &solver, const int vid,
//                    const vector<float> &distances, const float &max_dist);
vec3f polar_basis(const geodesic_solver &solver, const vector<vec3f> &positions,
                  const vector<vec3f> &normals, int vid);
Eigen::VectorXd wrapper(const vector<float> &f);

vector<float> exact_geodesic_distance(const vector<vec3i> &triangles,
                                      const vector<vec3f> &positions,
                                      const int &source);
std::tuple<Eigen::SparseMatrix<double, 1>, Eigen::SparseMatrix<double, 1>>
init_discrete_diff_op(
    const geodesic_solver &solver, const vector<vector<float>> &angles,
    const vector<float> &total_angles, const vector<vec3i> &triangles,
    const vector<vec3f> &positions, const vector<vec3i> &adjacencies,
    const vector<vector<int>> &v2t, const vector<vec3f> &normals);

void trace_in_triangles(const vector<vec3f> &positions,
                        const vector<vec3i> &triangles, const vec3f &dir,
                        const vec3f &bary, const int pid, vec3f &sample_pos,
                        vec3f &sample_bary);

int next_tid(const geodesic_solver &solver, const vector<vector<float>> &angles,
             const vector<vec3f> &positions, const vector<vector<int>> &v2t,
             const vector<vec3i> &triangles, const vector<vec3f> &normals,
             const int from, const vec3f &v);
pair<vector<mesh_point>, float> straightest_geodesic(
    const geodesic_solver &solver, const vector<vec3i> &triangles,
    const vector<vec3f> &positions, const vector<vec3f> &normals,
    const vector<vec3i> &adjacencies,
    const vector<vector<int>> &v2p_adjacencies,
    const vector<vector<float>> &angles, const vector<float> &total_angles,
    const mesh_point &from, const vec3f &v, const float &l,
    const int max_crossed_tri = INT_MAX);
void parallel_transp(const geodesic_solver &solver,
                     const vector<vector<float>> &angles,
                     const vector<float> &total_angles,
                     const vector<vec3i> &triangles,
                     const vector<vec3f> &positions,
                     const vector<vec3i> &adjacencies, vec3f &v,
                     const vector<vec3f> &normals, const int &from,
                     const int &to, const int &mode);

vector<vec3f> compute_grad(const geodesic_solver &solver,
                           const vector<vec3i> &triangles,
                           const vector<vec3f> &positions,
                           const vector<vec3f> &normals,
                           const Eigen::SparseMatrix<double, 1> &G,
                           const Eigen::VectorXd &f, bool normalized = false);

vector<vec3f> compute_grad(const geodesic_solver &solver,
                           const vector<vec3i> &triangles,
                           const vector<vec3f> &positions,
                           const vector<vec3f> &normals,
                           const Eigen::SparseMatrix<double, 1> &G,
                           const vector<float> &f, bool normalized = false);

vector<float> compute_laplacian(const discr_diff_op &operators,
                                const vector<float> &distances);

inline vec3f tri_bary_coords(const vector<vec3i> &triangles,
                             const vector<vec3f> &positions, const int pid,
                             const vec3f &p) {
  auto px = positions[triangles[pid].x];
  auto py = positions[triangles[pid].y];
  auto pz = positions[triangles[pid].z];
  return tri_bary_coords(px, py, pz, p);
}
inline vec3f tid_normal(const vector<vec3i> &triangles,
                        const vector<vec3f> &positions, const int tid) {
  auto p0 = positions[triangles[tid].x];
  auto p1 = positions[triangles[tid].y];
  auto p2 = positions[triangles[tid].z];

  return normalize(cross(p1 - p0, p2 - p0));
}

inline int node_is_adjacent(const geodesic_solver &solver, int vid, int node) {
  auto nbr = solver.graph[vid];
  for (auto i = 0; i < nbr.size(); ++i) {
    if (nbr[i].node == node) {
      return i;
    }
  }
  return -1;
}
template <typename T> inline int find(const T &vec, int x) {
  for (int i = 0; i < size(vec); i++)
    if (vec[i] == x)
      return i;
  return -1;
}