//
//  karcher.cpp
//  glfw
//
//  Created by Claudio Mancinelli on 21/07/2020.
//

#include "geometry.h"

#include <src/logging.h>
using namespace logging;
using namespace yocto;

#include <deque>

inline vec3f get_bary(const vec2f &uv) {
  return vec3f{1 - uv.x - uv.y, uv.x, uv.y};
}

// utility (geometry)
// project a vector v onto a plane having n as normal
vec3f project_vec(const vec3f &v, const vec3f &n) {
  auto proj = n * dot(v, n);

  return v - proj;
}
vec3f rot_vect(const vec3f &p, const vec3f &axis, const float angle) {
  auto M = rotation_frame(axis, angle);
  auto v = p;
  return transform_vector(M, v);
}
vec2f rot_vect(const vec2f &p, const float theta) {
  auto M = mat2f{{yocto::cos(theta), -yocto::sin(theta)},
                 {yocto::sin(theta), yocto::cos(theta)}};
  auto v = M * p;
  return v;
}

vec3f tid_centroid(const vector<vec3i> &triangles,
                   const vector<vec3f> &positions, const int tid) {
  vec3f p0 = positions[triangles[tid].x];
  vec3f p1 = positions[triangles[tid].y];
  vec3f p2 = positions[triangles[tid].z];

  return (p0 + p1 + p2) / 3.0;
}

mesh_point make_point_from_vert(const vector<vec3i> &triangles,
                                const vector<vector<int>> &v2t, const int vid) {
  auto tid = v2t[vid][0];
  auto k = find(triangles[tid], vid);
  auto bary = zero3f;
  bary[k] = 1.0;
  return {tid, {bary.y, bary.z}};
}
int forced_vert_from_point(const vector<vec3i> &triangles,
                           const mesh_point &p) {
  auto bary = vector<pair<float, int>>{
      {1 - p.uv.x - p.uv.y, 0}, {p.uv.x, 1}, {p.uv.y, 2}};

  sort(bary.begin(), bary.end());
  return triangles[p.face][bary.back().second];
}

void clean_bary(mesh_point &sample) {

  auto bary = sample.uv;
  auto coords = vector<pair<float, int>>{
      {1 - bary.x - bary.y, 0}, {bary.x, 1}, {bary.y, 2}};
  sort(coords.begin(), coords.end());
  vec3f bary3d = get_bary(bary);
  if (coords[0].first < 0 && coords[1].first > 0) {
    bary3d[coords[0].second] = 0;
    bary3d[coords[2].second] =
        1 - bary3d[coords[1].second] - bary3d[coords[0].second];
  } else if (coords[0].first < 0) {
    bary3d[coords[0].second] = 0;
    bary3d[coords[1].second] = 0;
    bary3d[coords[2].second] = 1;
  }

  sample = {sample.face, {bary3d.y, bary3d.z}};
}
float nbr_avg_edge_length(const geodesic_solver &G, int vid) {
  auto nbr = G.graph[vid];
  auto avg = 0.f;
  auto s = (int)nbr.size();
  for (int i = 0; i < s; ++i) {
    avg += nbr[i].length;
  }

  return (s != 0) ? avg / s : avg;
}

vec3f polar_basis(const geodesic_solver &solver, const vector<vec3f> &positions,
                  const vector<vec3f> &normals, int vid) {
  int vid0 = solver.graph[vid][0].node;
  vec3f v = positions[vid0] - positions[vid];
  vec3f e = project_vec(v, normals[vid]);
  return normalize(e);
}

vec3f polar_basis(const vector<vec3i> &triangles,
                  const vector<vec3f> &positions, int tid) {
  auto c = tid_centroid(triangles, positions, tid);
  vec3f v = positions[triangles[tid].x];
  return normalize(v - c);
}

// Compute polar coordinates
float angle_in_tangent_space(const geodesic_solver &solver,
                             const vector<vec3f> &positions, const vec3f &v,
                             const int vid, const vec3f &n) {
  float teta;
  int vid0 = solver.graph[vid][0].node;
  vec3f e0 = positions[vid0] - positions[vid];
  vec3f e = normalize(project_vec(e0, n));

  teta = angle(v, e);
  if (dot(cross(e, v), n) < 0)
    teta *= -1;

  return teta;
}
mesh_point point_from_vert(const vector<vec3i> &triangles,
                           const vector<vector<int>> &v2t, const int vid) {
  auto tid = v2t[vid][0];
  auto k = find(triangles[tid], vid);
  auto bary = zero3f;
  bary[k] = 1;
  return {tid, {bary.y, bary.z}};
}

vec3f tri_bary_coords(const vec3f &v0, const vec3f &v1, const vec3f &v2,
                      const vec3f &p) {
  vec3f wgts = vec3f{0.0, 0.0, 0.0};
  vec3f u = v1 - v0, v = v2 - v0, w = p - v0;
  float d00 = dot(u, u), d01 = dot(u, v), d11 = dot(v, v), d20 = dot(w, u),
        d21 = dot(w, v), d = d00 * d11 - d01 * d01;

  if (d == 0)
    return zero3f;

  wgts[2] = (d00 * d21 - d01 * d20) / d;
  assert(!isnan(wgts[2]));
  wgts[1] = (d11 * d20 - d01 * d21) / d;
  assert(!isnan(wgts[1]));
  wgts[0] = 1.0 - wgts[1] - wgts[2];
  assert(!isnan(wgts[0]));

  return wgts;
}

// returns the index of the triangles in the star of from which v is pointing
// to
int next_tid(const geodesic_solver &solver, const vector<vector<float>> &angles,
             const vector<vec3f> &positions, const vector<vector<int>> &v2t,
             const vector<vec3i> &triangles, const vector<vec3f> &normals,
             const int from, const vec3f &v) {
  auto teta = angle_in_tangent_space(solver, positions, v, from, normals[from]);
  if (teta < 0)
    teta += 2 * M_PI;
  auto nbr = angles[from];
  int s = (int)nbr.size();
  if (teta == 0)
    return v2t[from][0];
  for (int i = 0; i < s; ++i) {
    if (nbr[i] < teta)
      continue;

    if (i % 2 == 0) {
      return v2t[from][(i - 2) / 2];
    } else {
      return v2t[from][(i - 1) / 2];
    }
  }
  return v2t[from].back();
}
vec3f trace_segment_vert(const vector<vec3f> &verts, const vec3f n,
                         const vec3f bary, const vec3f baryV, const vec3f &dir,
                         const int offset) {
  auto right = verts[(offset + 1) % 3] - verts[offset];
  auto left = verts[(offset + 2) % 3] - verts[offset];
  auto sample_bary = zero3f;
  if (dot(cross(right, dir), n) > 0) {
    if (dot(cross(dir, left), n) > 0) {
      auto factor = bary[offset] / baryV[offset];
      sample_bary[offset] = 0;
      sample_bary[(offset + 1) % 3] =
          bary[(offset + 1) % 3] - baryV[(offset + 1) % 3] * factor;
      sample_bary[(offset + 2) % 3] =
          bary[(offset + 2) % 3] - baryV[(offset + 2) % 3] * factor;
    } else
      sample_bary[(offset + 2) % 3] = 1;
  } else
    sample_bary[(offset + 1) % 3] = 1;

  return sample_bary;
}

vec3f trace_segment_edge(const vector<vec3f> &verts, const vec3f n,
                         const vec3f bary, const vec3f baryV, const vec3f &dir,
                         const int offset, const vec3f &sample_coords) {
  auto sample_bary = zero3f;
  auto right = verts[(offset + 1) % 3] - sample_coords;
  auto left = verts[offset] - sample_coords;
  auto front = verts[(offset + 2) % 3] - sample_coords;
  if (dot(cross(right, dir), n) > 0) {
    if (dot(cross(dir, front), n) > 0) {
      auto factor = bary[offset] / baryV[offset];
      sample_bary[offset] = 0;
      sample_bary[(offset + 1) % 3] =
          bary[(offset + 1) % 3] - baryV[(offset + 1) % 3] * factor;
      sample_bary[(offset + 2) % 3] =
          bary[(offset + 2) % 3] - baryV[(offset + 2) % 3] * factor;

    } else {
      if (dot(cross(dir, left), n) > 0) {
        auto factor = bary[(offset + 1) % 3] / baryV[(offset + 1) % 3];
        sample_bary[(offset + 1) % 3] = 0;
        sample_bary[offset] = bary[offset] - baryV[offset] * factor;
        sample_bary[(offset + 2) % 3] =
            bary[(offset + 2) % 3] - baryV[(offset + 2) % 3] * factor;

      } else {
        if (dot(left, dir) > 0) {
          sample_bary[(offset)] = 1;
        } else {
          sample_bary[(offset + 1) % 3] = 1;
        }
      }
    }
  } else {
    if (dot(right, dir) > 0) {
      sample_bary[(offset + 1) % 3] = 1;
    } else {
      sample_bary[(offset)] = 1;
    }
  }

  return sample_bary;
}
vec3f trace_segment_tri(const vector<vec3f> &verts, const vec3f n,
                        const vec3f bary, const vec3f baryV, const vec3f &dir,
                        const vec3f &sample_coords) {
  auto sample_bary = zero3f;
  vec3f w0 = verts[0] - sample_coords;
  vec3f w1 = verts[1] - sample_coords;
  vec3f w2 = verts[2] - sample_coords;
  if (dot(cross(w0, dir), n) > 0 && dot(cross(dir, w1), n) > 0) {
    sample_bary = vec3f{bary[0] - bary[2] * baryV[0] / baryV[2],
                        bary[1] - bary[2] * baryV[1] / baryV[2], 0};

  } else if (dot(cross(w1, dir), n) > 0 && dot(cross(dir, w2), n) > 0) {
    sample_bary = vec3f{0, bary[1] - bary[0] * baryV[1] / baryV[0],
                        bary[2] - bary[0] * baryV[2] / baryV[0]};
  } else {
    sample_bary = vec3f{bary[0] - bary[1] * baryV[0] / baryV[1], 0,
                        bary[2] - bary[1] * baryV[2] / baryV[1]};
  }

  return sample_bary;
}
// Identify the intersection of the polyline inside triangle pid
// tracing: https://cims.nyu.edu/gcl/papers/campen2016bms.pdf
void trace_in_triangles(const vector<vec3f> &positions,
                        const vector<vec3i> &triangles, const vec3f &dir,
                        const vec3f &bary, const int pid, vec3f &sample_pos,
                        vec3f &sample_bary) {
  vec3f baryM = zero3f, baryV = zero3f;
  vec3f v0 = positions[triangles[pid].x];
  vec3f v1 = positions[triangles[pid].y];
  vec3f v2 = positions[triangles[pid].z];
  vector<vec3f> verts = {v0, v1, v2};
  vec3f n = triangle_normal(v0, v1, v2);
  vec3f sample_coords = bary.x * v0 + bary.y * v1 + bary.z * v2;
  vec3f M = sample_coords + dir;

  baryM = tri_bary_coords(v0, v1, v2, M);
  for (int i = 0; i < 3; ++i) {
    baryV[i] = baryM[i] - bary[i];
  }
  auto [is_vertex, k_vert] = bary_is_vert(bary);
  auto [is_on_edge, k_edge] = bary_is_edge(bary);
  if (is_vertex) {
    sample_bary = trace_segment_vert(verts, n, bary, baryV, dir, k_vert);
    sample_pos = sample_bary.x * verts[0] + sample_bary.y * verts[1] +
                 sample_bary.z * verts[2];
  } else if (is_on_edge) {
    sample_bary =
        trace_segment_edge(verts, n, bary, baryV, dir, k_edge, sample_coords);
    sample_pos = sample_bary.x * verts[0] + sample_bary.y * verts[1] +
                 sample_bary.z * verts[2];
  } else {
    sample_bary = trace_segment_tri(verts, n, bary, baryV, dir, sample_coords);
    sample_pos = sample_bary.x * verts[0] + sample_bary.y * verts[1] +
                 sample_bary.z * verts[2];
  }
}

void parallel_transp(const geodesic_solver &solver,
                     const vector<vector<float>> &angles,
                     const vector<float> &total_angles,
                     const vector<vec3i> &triangles,
                     const vector<vec3f> &positions,
                     const vector<vec3i> &adjacencies, vec3f &v,
                     const vector<vec3f> &normals, const int &from,
                     const int &to, const int &mode) {
  switch (mode) {
  case V2V: {
    auto mag = length(v);
    float teta =
        angle_in_tangent_space(solver, positions, v, from, normals[from]);
    auto nbr_from = solver.graph[from];
    auto nbr_to = solver.graph[to];
    float phi_ij = -1;
    float phi_ji = -1;

    for (int i = 0; i < nbr_from.size(); ++i) {
      if (nbr_from[i].node == to) {
        phi_ij = angles[from][i];

        break;
      }
    }

    for (int j = 0; j < nbr_to.size(); ++j) {
      if (nbr_to[j].node == from) {
        phi_ji = angles[to][j];
        break;
      }
    }
    assert(phi_ij != -1);
    assert(phi_ji != -1);

    vec3f e0 = polar_basis(solver, positions, normals, to);
    float rotation = teta + phi_ji + M_PI - phi_ij;

    v = rot_vect(e0, normals[to], rotation);
    v *= mag;

  } break;

  case V2T: {
    float teta =
        angle_in_tangent_space(solver, positions, v, from, normals[from]);
    vec3i tri = triangles[to];
    vec3f p0 = positions[tri.x];
    vec3f p1 = positions[tri.y];
    vec3f p2 = positions[tri.z];
    vec3f normal = triangle_normal(p0, p1, p2);
    vec3f centroid = (p0 + p1 + p2) / 3.0;
    vec3f e = normalize(p0 - centroid);

    vec3f coords = positions[from] - centroid;
    float phi_ji = angle(e, coords);
    if (dot(cross(e, coords), normal) < 0)
      phi_ji = 2 * M_PI - phi_ji;
    int offset = find(tri, from);
    assert(offset != -1);
    int vid1 = tri[(offset + 1) % 3];
    int vid2 = tri[(offset + 2) % 3];
    float factor = 2 * M_PI / total_angles[from];
    auto nbr_from = solver.graph[from];
    float phi_ij = -1;
    coords *= -1;
    if (nbr_from[0].node == vid2) {
      vec3f edge = positions[vid2] - positions[from];
      float curr_angle = angle(edge, coords);
      curr_angle *= factor;
      curr_angle = 2 * M_PI - curr_angle;
      phi_ij = curr_angle;
    } else {
      for (int i = 0; i < nbr_from.size(); ++i) {
        if (nbr_from[i].node == vid1) {
          phi_ij = angles[from][i];
          break;
        }
      }

      vec3f edge = positions[vid1] - positions[from];
      float curr_angle = angle(edge, coords);
      curr_angle *= factor;
      phi_ij += curr_angle;
    }

    float rot = teta + phi_ji + M_PI - phi_ij;

    e *= length(v);
    v = rot_vect(e, normal, rot);

  }

  break;

  case T2V: {
    vec3i tri = triangles[from];
    vec3f p0 = positions[tri.x];
    vec3f p1 = positions[tri.y];
    vec3f p2 = positions[tri.z];
    vec3f n = triangle_normal(p0, p1, p2);
    vec3f centroid = (p0 + p1 + p2) / 3.0;
    vec3f e = normalize(p0 - centroid);
    float teta = angle(e, v);

    if (dot(cross(e, v), n) < 0)
      teta = 2 * M_PI - teta;
    int offset = find(tri, to);
    assert(offset != -1);
    int vid1 = tri[(offset + 1) % 3];

    vec3f vert = positions[tri[offset]];
    vec3f v1 = positions[vid1] - vert;

    vec3f coords = vert - centroid;
    float phi_ij = angle(e, coords);
    if (dot(cross(e, coords), n) < 0)
      phi_ij = 2 * M_PI - phi_ij;

    coords *= -1;
    float phi_ji = angle(v1, coords);
    float factor = 2 * M_PI / total_angles[to];
    phi_ji *= factor;
    auto nbr = solver.graph[to];
    for (int i = 0; i < nbr.size(); ++i) {
      if (nbr[i].node == vid1) {
        float phi = angles[to][i];
        phi_ji += phi;
        break;
      }
    }

    float rot = teta + phi_ji + M_PI - phi_ij;
    vec3f e0 = polar_basis(solver, positions, normals, to);
    e0 *= length(v);
    v = rot_vect(e0, normals[to], rot);

  } break;

  case T2T: {
    auto flat_from = init_flat_triangle(positions, triangles[from]);
    auto k = find(adjacencies[from], to);
    assert(k != -1);
    auto flat_to = unfold_face(triangles, positions, flat_from, from, to);
    auto bary = vec2f{0.333, 0.333};
    auto c0 =
        interpolate_triangle(flat_from[0], flat_from[1], flat_from[2], bary);
    auto c1 = interpolate_triangle(flat_to[0], flat_to[1], flat_to[2], bary);
    auto e0 = flat_from[0] - c0;
    auto e1 = flat_to[0] - c1;

    auto w = c1 - c0;
    auto phi_ij = angle(e0, w);
    if (cross(e0, w) < 0)
      phi_ij = 2 * M_PI - phi_ij;
    w *= -1;
    auto phi_ji = angle(e1, w);
    if (cross(e1, w) < 0)
      phi_ji = 2 * M_PI - phi_ji;

    auto n = tid_normal(triangles, positions, from);
    auto e = polar_basis(triangles, positions, from);
    float teta = angle(e, v);
    if (dot(cross(e, v), n) < 0)
      teta = 2 * M_PI - teta;

    auto e_to = polar_basis(triangles, positions, to);
    auto n_to = tid_normal(triangles, positions, to);
    float rot = teta + phi_ji + M_PI - phi_ij;
    e_to *= length(v);
    v = rot_vect(e_to, n_to, rot);

  }

  break;
  }
}
// VTP
vector<float> exact_geodesic_distance(const vector<vec3i> &triangles,
                                      const vector<vec3f> &positions,
                                      const int &source) {
  int V = (int)positions.size();
  int F = (int)triangles.size();
  vector<double> points(3 * V);
  vector<uint> faces(3 * F);
  vector<float> f(V);

  for (int i = 0; i < V; ++i) {
    points[3 * i] = positions[i].x;
    points[3 * i + 1] = positions[i].y;
    points[3 * i + 2] = positions[i].z;
  }
  for (int i = 0; i < F; ++i) {
    faces[3 * i] = triangles[i].x;
    faces[3 * i + 1] = triangles[i].y;
    faces[3 * i + 2] = triangles[i].z;
  }
  geodesic_VTP::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);
  geodesic_VTP::GeodesicAlgorithmExact algorithm(&mesh);
  algorithm.propagate(source);
  vector<geodesic_VTP::Vertex> verts = mesh.vertices();
  for (int j = 0; j < V; ++j) {
    geodesic_VTP::Vertex v = verts[j];
    float value = (float)v.geodesic_distance();
    f[j] = value;
  }

  return f;
}

// utility (gradient matrix)
Eigen::VectorXd wrapper(const vector<float> &f) {
  Eigen::VectorXd F(f.size());
  for (int i = 0; i < f.size(); ++i) {
    F(i) = f[i];
  }
  return F;
}

Eigen::MatrixXd rhs(int s) {
  Eigen::MatrixXd E(s, s + 1);
  Eigen::MatrixXd X = Eigen::MatrixXd::Constant(s, 1, -1);
  E.topLeftCorner(s, 1) = X;
  Eigen::MatrixXd I(s, s);
  I.setIdentity();
  E.topRightCorner(s, s) = I;
  return E;
}
void fill_riemannian_gradient_entries(
    vector<Eigen::Triplet<double>> &entries,
    const vector<vector<pair<int, float>>> &ring, const Eigen::VectorXd &c,
    const Eigen::VectorXd &a0, const Eigen::VectorXd &a1, const int n) {
  int vid = ring[0][0].first;
  int s = (int)ring.size();
  double c0_squared = pow(c[0], 2);
  double c1_squared = pow(c[1], 2);
  Eigen::Matrix2d g_inv;
  double det = 1 + c0_squared + c1_squared;
  g_inv << 1 + c1_squared, -c[0] * c[1], -c[0] * c[1], 1 + c0_squared;
  g_inv /= det;
  typedef Eigen::Triplet<double> T;
  for (int i = 0; i < s; ++i) {
    for (auto j = 0; j < ring[i].size(); ++j) {
      int entry = ring[i][j].first;
      auto w = ring[i][j].second;
      entries.push_back(
          T(vid, entry, w * (g_inv(0, 0) * a0(i) + g_inv(0, 1) * a1(i))));
      entries.push_back(
          T(n + vid, entry, w * (g_inv(1, 0) * a0(i) + g_inv(1, 1) * a1(i))));
    }
  }
}

void laplacian_entries(vector<Eigen::Triplet<double>> &entries,
                       const Eigen::Matrix2d &g, const Eigen::Matrix2d &g_inv,
                       const vector<vector<pair<int, float>>> &ring,
                       const Eigen::VectorXd &c, const Eigen::MatrixXd &a) {
  int vid = ring[0][0].first;
  typedef Eigen::Triplet<double> T;
  Eigen::VectorXd b;
  double c0_squared = pow(c[0], 2);
  double c1_squared = pow(c[1], 2);

  double det = 1 + c0_squared + c1_squared;

  double g11u = 2 * c[0] * c[2], g12u = c[0] * c[3] + c[1] * c[2],
         g22_u = 2 * c[1] * c[3], g11_v = 2 * c[0] * c[3],
         g12_v = c[0] * c[4] + c[1] * c[3], g22_v = 2 * c[1] * c[4],
         g_u = g11u + g22_u, g_v = g11_v + g22_v,
         g_invu11 = (det * g22_u - g(1, 1) * g_u) / pow(det, 2),
         g_invu12 = -(det * g12u - g(0, 1) * g_u) / pow(det, 2),
         g_invv21 = -(det * g12_v - g(0, 1) * g_v) / pow(det, 2),
         g_invv22 = (det * g11_v - g(0, 0) * g_v) / pow(det, 2);

  double coeff0 = (g_u * g_inv(0, 0) + g_v * g_inv(1, 0)) / (2 * det) +
                  g_invu11 + g_invv21,
         coeff1 = (g_u * g_inv(0, 1) + g_v * g_inv(1, 1)) / (2 * det) +
                  g_invu12 + g_invv22,
         coeff2 = g_inv(0, 0), coeff3 = 2 * g_inv(0, 1), coeff4 = g_inv(1, 1);

  b = coeff0 * a.row(0) + coeff1 * a.row(1) + coeff2 * a.row(2) +
      coeff3 * a.row(3) + coeff4 * a.row(4);
  for (int i = 0; i < b.size(); ++i) {
    for (auto j = 0; j < ring[i].size(); ++j) {
      int entry = ring[i][j].first;

      entries.push_back(T(vid, entry, b(i) * ring[i][j].second));
    }
  }
}

std::tuple<Eigen::SparseMatrix<double, 1>, Eigen::SparseMatrix<double, 1>>
init_discrete_diff_op(
    const geodesic_solver &solver, const vector<vector<float>> &angles,
    const vector<float> &total_angles, const vector<vec3i> &triangles,
    const vector<vec3f> &positions, const vector<vec3i> &adjacencies,
    const vector<vector<int>> &v2t, const vector<vec3f> &normals) {

  typedef Eigen::Triplet<double> T;
  vector<T> L_entries;
  vector<T> G_entries;
  Eigen::SparseMatrix<double, 1> Grad;
  Eigen::SparseMatrix<double, 1> Lap;
  int V = (int)positions.size();
  for (int i = 0; i < V; ++i) {

    auto nbr = solver.graph[i];
    vec3f vert = positions[i];
    vec3f n = normals[i];
    int s = (int)nbr.size();
    Eigen::MatrixXd Q(s, 5);
    Eigen::VectorXd h(s);
    auto teta = 0.f;
    auto pos = zero2f;
    auto d = 0.f;
    auto coords = zero3f;
    vector<vector<pair<int, float>>> a_map(s + 1);
    a_map[0].push_back({i, 1});
    auto len = nbr_avg_edge_length(solver, i);
    auto e = polar_basis(solver, positions, normals, i);
    auto start = point_from_vert(triangles, v2t, i);
    for (int j = 0; j < s; ++j) {
      int curr = nbr[j].node;
      if (j % 2 == 1) {

        auto teta_prev = angles[i][j - 1];
        auto teta_next = (j == s - 1) ? 2 * pif : angles[i][j + 1];
        auto curr_teta = (teta_prev + teta_next) / 2;

        auto dir = rot_vect(e, n, curr_teta);
        auto path = straightest_geodesic(solver, triangles, positions, normals,
                                         adjacencies, v2t, angles, total_angles,
                                         start, dir, len, 2);
        auto sample = path.first.back();
        auto point = eval_position(triangles, positions, sample);
        auto curr_len = path.second;
        pos = vec2f{curr_len * yocto::cos(curr_teta),
                    curr_len * yocto::sin(curr_teta)};
        coords = point - vert;
        auto bary_end = get_bary(sample.uv);
        for (auto h = 0; h < 3; ++h) {
          a_map[j + 1].push_back({triangles[sample.face][h], bary_end[h]});
        }

      } else {

        a_map[j + 1].push_back({nbr[j].node, 1});
        teta = angles[i][j];
        d = nbr[j].length;
        pos = vec2f{d * std::cos(teta), d * std::sin(teta)};
        coords = positions[curr] - vert;
      }

      Q(j, 0) = pos[0];
      Q(j, 1) = pos[1];
      Q(j, 2) = pow(pos[0], 2) / 2;
      Q(j, 3) = pos[0] * pos[1];
      Q(j, 4) = pow(pos[1], 2) / 2;

      h(j) = dot(coords, n);
    }
    Eigen::MatrixXd Qt = Eigen::Transpose<Eigen::MatrixXd>(Q);

    Eigen::MatrixXd A = Qt * Q;
    Eigen::MatrixXd E = rhs(s);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(A);
    Eigen::VectorXd c(5);
    Eigen::MatrixXd a(5, s + 1);

    if (dec.isInvertible()) {
      Eigen::MatrixXd inv = A.inverse();
      c = inv * Qt * h;
      a = inv * Qt * E;

    } else {
      Eigen::MatrixXd Rhsc = Qt * h;
      Eigen::MatrixXd Rhsa = Qt * E;
      c = dec.solve(Rhsc);
      a = dec.solve(Rhsa);
    }
    Eigen::Matrix2d g;
    Eigen::Matrix2d g_inv;

    double c0_squared = pow(c[0], 2);
    double c1_squared = pow(c[1], 2);
    g << 1 + c0_squared, c[0] * c[1], c[0] * c[1], 1 + c1_squared;

    double det = 1 + c0_squared + c1_squared;
    g_inv << 1 + c1_squared, -c[0] * c[1], -c[0] * c[1], 1 + c0_squared;
    g_inv /= det;

    laplacian_entries(L_entries, g, g_inv, a_map, c, a);

    fill_riemannian_gradient_entries(G_entries, a_map, c, a.row(0), a.row(1),
                                     V);
  }

  Lap.resize(V, V);
  Grad.resize(2 * V, V);
  Lap.setFromTriplets(L_entries.begin(), L_entries.end());
  Grad.setFromTriplets(G_entries.begin(), G_entries.end());

  return {Grad, Lap};
}

vec3f polar_to_cartesian(const geodesic_solver &solver,
                         const vector<vec3f> &positions,
                         const vector<vec3f> &normals, const double x,
                         const double y, const int vid) {
  vec3f g = zero3f;
  vec2f sol = vec2f{(float)x, (float)y};
  double phi = yocto::atan2(y, x);

  float mag = length(sol);
  vec3f e = polar_basis(solver, positions, normals, vid);
  g = rot_vect(e, normals[vid], phi);
  g *= mag;

  return g;
}
vector<vec3f> compute_grad(const geodesic_solver &solver,
                           const vector<vec3i> &triangles,
                           const vector<vec3f> &positions,
                           const vector<vec3f> &normals,
                           const Eigen::SparseMatrix<double, 1> &G,
                           const Eigen::VectorXd &f, bool normalized) {
  auto V = (int)positions.size();

  auto F = (int)triangles.size();

  if (G.rows() == 2 * V) {
    vector<vec3f> g(V);

    Eigen::VectorXd Grad = G * f;

    for (int i = 0; i < V; ++i) {
      g[i] = polar_to_cartesian(solver, positions, normals, Grad(i),
                                Grad(V + i), i);
      if (normalized)
        g[i] = normalize(g[i]);
    }
    return g;
  } else {
    vector<vec3f> g(F);
    Eigen::VectorXd Grad = G * f;

    for (int i = 0; i < F; ++i) {
      g[i].x = Grad(3 * i);
      g[i].y = Grad(3 * i + 1);
      g[i].z = Grad(3 * i + 2);

      if (normalized)
        g[i] = normalize(g[i]);
    }
    return g;
  }
}
vector<vec3f> compute_grad(const geodesic_solver &solver,
                           const vector<vec3i> &triangles,
                           const vector<vec3f> &positions,
                           const vector<vec3f> &normals,
                           const Eigen::SparseMatrix<double, 1> &G,
                           const vector<float> &f, bool normalized) {
  auto F = wrapper(f);
  return compute_grad(solver, triangles, positions, normals, G, F, normalized);
}
vector<float> compute_laplacian(const discr_diff_op &operators,
                                const vector<float> &distances) {

  auto DIST = wrapper(distances);
  vector<float> laplacian(distances.size());
  Eigen::VectorXd Lap = operators.Lap * DIST;
  for (auto i = 0; i < Lap.size(); ++i) {
    laplacian[i] = Lap(i);
  }
  return laplacian;
}

pair<int, vec3f> handle_vert(const geodesic_solver &solver,
                             const vector<vec3i> &triangles,
                             const vector<vec3f> &positions,
                             const vector<vec3f> &normals,
                             const vector<vec3i> &adjacencies,
                             const vector<vector<int>> &v2p_adjacencies,
                             const vector<vector<float>> &angles,
                             const vector<float> &total_angles, const int vid,
                             const int tid, const vec3f dir) {
  auto v = dir;
  parallel_transp(solver, angles, total_angles, triangles, positions,
                  adjacencies, v, normals, tid, vid, T2V);
  auto next = next_tid(solver, angles, positions, v2p_adjacencies, triangles,
                       normals, vid, v);

  parallel_transp(solver, angles, total_angles, triangles, positions,
                  adjacencies, v, normals, vid, next, V2T);
  return std::make_pair(next, v);
}

pair<vector<mesh_point>, float> straightest_geodesic(
    const geodesic_solver &solver, const vector<vec3i> &triangles,
    const vector<vec3f> &positions, const vector<vec3f> &normals,
    const vector<vec3i> &adjacencies,
    const vector<vector<int>> &v2p_adjacencies,
    const vector<vector<float>> &angles, const vector<float> &total_angles,
    const mesh_point &from, const vec3f &v, const float &l,
    const int max_crossed_tri) {
  auto vid = -1, tid = -1, next_tri = -1;
  auto dir = v;
  float len = 0.0;
  auto next_bary = zero3f, next = zero3f;
  auto prev = eval_position(triangles, positions, from);
  auto samples = vector<mesh_point>{from};
  auto bary = get_bary(from.uv);
  auto [is_vert, kv] = bary_is_vert(bary);
  auto [is_on_edge, ke] = bary_is_edge(bary);
  auto crossed_tri = 0;
  if (is_vert) {
    vid = triangles[from.face][kv];
    tid = next_tid(solver, angles, positions, v2p_adjacencies, triangles,
                   normals, vid, v);
    kv = find(triangles[tid], vid);
    bary = zero3f;
    bary[kv] = 1;
    parallel_transp(solver, angles, total_angles, triangles, positions,
                    adjacencies, dir, normals, vid, tid, V2T);
  } else if (is_on_edge) {
    auto p0 = triangles[from.face][ke];
    auto p1 = triangles[from.face][(ke + 1) % 3];
    auto p2 = triangles[from.face][(ke + 2) % 3];
    auto n = triangle_normal(positions[p0], positions[p1], positions[p2]);
    auto edge = normalize(positions[p1] - positions[p0]);
    if (dot(cross(edge, v), n) > 0)
      tid = from.face;
    else {
      tid = adjacencies[from.face][ke];
      bary = tri_bary_coords(triangles, positions, tid, prev);
      parallel_transp(solver, angles, total_angles, triangles, positions,
                      adjacencies, dir, normals, from.face, tid, T2T);
    }

  } else
    tid = from.face;

  while (len < l && crossed_tri < max_crossed_tri) {
    trace_in_triangles(positions, triangles, dir, bary, tid, next, next_bary);
    samples.push_back({tid, vec2f{next_bary.y, next_bary.z}});
    len += length(next - prev);
    ++crossed_tri;
    if (len < l) {
      prev = next;
      auto [V, k_v] = bary_is_vert(next_bary);
      auto [E, k_e] = bary_is_edge(next_bary);
      if (V) {
        vid = triangles[tid][k_v];
        auto out =
            handle_vert(solver, triangles, positions, normals, adjacencies,
                        v2p_adjacencies, angles, total_angles, vid, tid, dir);
        tid = out.first;
        dir = out.second;
        k_v = find(triangles[tid], vid);
        bary = zero3f;
        bary[k_v] = 1;
      } else if (E) {
        next_tri = adjacencies[tid][k_e];
        auto p0 = triangles[tid][k_e];
        auto offset0 = find(triangles[next_tri], p0);
        auto offset1 = (offset0 + 2) % 3;
        bary = zero3f;
        bary[offset0] = next_bary[k_e];
        bary[offset1] = next_bary[(k_e + 1) % 3];

        parallel_transp(solver, angles, total_angles, triangles, positions,
                        adjacencies, dir, normals, tid, next_tri, T2T);
        tid = next_tri;
      } else
        assert(false);
    }
  }
  if (len > l) {
    auto factor = (len - l);
    auto w = normalize(prev - next);
    w *= factor;
    w += next;
    bary = tri_bary_coords(triangles, positions, tid, w);
    samples.pop_back();
    samples.push_back({tid, vec2f{bary.y, bary.z}});
    len = l;
  } else {
    samples.push_back({tid, {bary.y, bary.z}});
  }

  return {samples, len};
}
std::tuple<vector<int>, mesh_point, mesh_point>
handle_short_strips(const vector<vec3i> &triangles,
                    const vector<vec3f> &positions, const vector<int> &strip,
                    const mesh_point &start, const mesh_point &end) {
  if (strip.size() == 1) {
    return {strip, start, end};
  } else if (strip.size() == 2) {
    auto [inside, b2f] =
        point_in_triangle(triangles, positions, start.face,
                          eval_position(triangles, positions, end));
    if (inside) {
      auto new_end = mesh_point{start.face, b2f};
      return {{start.face}, start, new_end};
    }
    std::tie(inside, b2f) =
        point_in_triangle(triangles, positions, end.face,
                          eval_position(triangles, positions, start));

    if (inside) {
      auto new_start = mesh_point{end.face, b2f};
      return {{end.face}, new_start, end};
    }

    return {strip, start, end};
  }
  return {{-1}, {}, {}};
}
vec3f flip_bary_to_adjacent_tri(const vector<vec3i> &adjacencies,
                                const int tid0, const int tid1,
                                const vec3f &bary) {
  if (tid0 == tid1)
    return bary;
  auto new_bary = zero3f;
  auto k1 = find(adjacencies[tid1], tid0);
  auto k0 = find(adjacencies[tid0], tid1);
  if (k1 == -1) {
    std::cout << "Error, faces are not adjacent" << std::endl;
    return zero3f;
  }
  new_bary[k1] = bary[(k0 + 1) % 3];
  new_bary[(k1 + 1) % 3] = bary[k0];
  new_bary[(k1 + 2) % 3] = bary[(k0 + 2) % 3];

  return new_bary;
}
