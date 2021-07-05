#include "cut_locus.h"

#include "strip.h"
#include <deque>
#include <mutex>

// bool is_edge(const shape_data &mesh, mesh_point &start, mesh_point &end) {
//   auto [start_is_vert, k_start] = point_is_vert(start);
//   auto [end_is_vert, k_end] = point_is_vert(end);
//   if (!(start_is_vert && end_is_vert))
//     return false;
//   auto start_vid = mesh.triangles[start.face][k_start];
//   auto end_vid = mesh.triangles[end.face][k_end];
//   auto entry = node_is_adjacent(mesh.solver, start_vid, end_vid);
//   if (entry % 2 != 0)
//     return false;
//   auto tid = mesh.v2t[start_vid][entry / 2];
//   auto ks = find(mesh.triangles[tid], start_vid);
//   auto ke = find(mesh.triangles[tid], end_vid);
//   auto bary = zero3f;
//   bary[ks] = 1;
//   start = {tid, {bary.y, bary.z}};
//   bary = zero3f;
//   bary[ke] = 1;
//   end = {tid, {bary.y, bary.z}};
//   return true;
// }

int field_maximum(const vector<float> &distances) {

  auto lambda = flt_min;
  auto max = -1;
  for (auto i = 0; i < distances.size(); ++i) {
    if (distances[i] > lambda) {
      lambda = distances[i];
      max = i;
    }
  }

  return max;
}
geodesic_path compute_geodesic_path(const shape_data &mesh,
                                    const shape_topology &topology,
                                    const dual_geodesic_solver &dual_solver,
                                    const geodesic_solver &solver,
                                    const mesh_point &start,
                                    const mesh_point &end) {

  auto path = geodesic_path{};
  if (start.face == end.face) {
    path.start = start;
    path.end = end;
    path.lerps = {};
    path.strip = {start.face};
    return path;
  }
  auto parents = vector<int>{};
  auto strip = compute_strip_tlv(mesh, dual_solver, end.face, start.face);

  path = shortest_path(mesh.triangles, mesh.positions, topology.adjacencies,
                       start, end, strip);
  return path;
}
std::tuple<vector<pair<float, int>>, int, float>
compute_laplacian_with_extrema(const discr_diff_op &operators,
                               const vector<float> &distances) {
  auto lap = compute_laplacian(operators, distances);
  auto result = vector<pair<float, int>>(lap.size());
  auto min = -1;
  auto lambda = flt_max;
  auto max = flt_min;

  for (auto i = 0; i < lap.size(); ++i) {
    result[i] = {lap[i], i};
    if (lap[i] < lambda) {
      min = i;
      lambda = lap[i];
    }
    if (max < yocto::abs(lap[i]))
      max = yocto::abs(lap[i]);
  }

  return {result, min, max};
}

float deviation_of_gradient(const shape_data &mesh,
                            const shape_topology &topology,
                            const geodesic_solver &solver,
                            const vector<vec3f> &gradients,
                            const int curr_node) {

  auto ref = gradients[curr_node];
  auto n = mesh.normals[curr_node];
  auto teta = angle_in_tangent_space(solver, mesh.positions, ref, curr_node, n);
  auto nbr = solver.graph[curr_node];
  auto lambda = flt_min;
  auto phi = 0.f;
  if (teta < 0)
    teta += 2 * pif;

  if (teta <= pif) {
    for (auto i = 0; i < nbr.size(); ++i) {
      if (i % 2 == 1)
        continue;
      auto v = gradients[nbr[i].node];
      if (length(v) <= 1e-1)
        continue;
      parallel_transport(mesh, topology, solver, v, nbr[i].node, curr_node,
                         V2V);
      phi = angle_in_tangent_space(solver, mesh.positions, v, curr_node, n);
      if (phi < 0)
        phi += 2 * pif;
      if (topology.angles[curr_node][i] <= teta ||
          topology.angles[curr_node][i] >= pif + teta) {
        if (phi <= teta || phi >= pif + teta)
          phi = -angle(v, ref);
        else
          phi = angle(v, ref);
      } else {
        if (phi > teta && phi < pif + teta)
          phi = -angle(v, ref);
        else
          phi = angle(v, ref);
      }
      if (lambda < phi)
        lambda = phi;
    }
  } else {
    for (auto i = 0; i < nbr.size(); ++i) {
      if (i % 2 == 1)
        continue;
      auto v = gradients[nbr[i].node];
      if (length(v) <= 1e-1)
        continue;
      parallel_transport(mesh, topology, solver, v, nbr[i].node, curr_node,
                         V2V);
      phi = angle_in_tangent_space(solver, mesh.positions, v, curr_node, n);
      if (phi < 0)
        phi += 2 * pif;
      if (topology.angles[curr_node][i] >= teta ||
          topology.angles[curr_node][i] <= teta - pif) {
        if (phi >= teta || phi <= teta - pif)
          phi = -angle(v, ref);
        else
          phi = angle(v, ref);
      } else {
        if (phi > teta - pif && phi < teta)
          phi = -angle(v, ref);
        else
          phi = angle(v, ref);
      }
      if (lambda < phi)
        lambda = phi;
    }
  }
  return lambda;
}

void update_angle_deviation(const shape_data &mesh,
                            const shape_topology &topology,
                            const geodesic_solver &solver,
                            const vector<vec3f> &gradients,
                            vector<node> &graph) {
  time_function();
  graph.resize(mesh.positions.size(), node{});
  for (auto i = 0; i < mesh.positions.size(); ++i) {
    graph[i].angle_deviation =
        deviation_of_gradient(mesh, topology, solver, gradients, i);
  }
}

vec2i compute_parent_along_cotree(const shape_data &mesh,
                                  const shape_topology &topology,
                                  const vector<node> &cotree,
                                  const vector<node> &graph, const vec2i &gen) {

  auto coparent0x = cotree[gen.x].neighbors[0];
  auto coparent1x = cotree[coparent0x].neighbors[0];
  auto coparent0y = cotree[gen.y].neighbors[0];
  auto coparent1y = cotree[coparent0y].neighbors[0];
  bool stopx = false;
  bool stopy = false;
  auto has_been_crossed = vector<bool>(mesh.triangles.size(), false);
  while (true) {
    if (!stopx) {
      auto k = find(topology.adjacencies[coparent0x], coparent1x);
      auto primal_x = mesh.triangles[coparent0x][k];
      auto primal_y = mesh.triangles[coparent0x][(k + 1) % 3];
      has_been_crossed[coparent0x] = true;
      has_been_crossed[coparent1x] = true;
      if (graph[primal_x].is_in_restricted && graph[primal_y].is_in_restricted)
        return vec2i{primal_x, primal_y};
    }

    if (!stopy) {
      auto k = find(topology.adjacencies[coparent0y], coparent1y);
      auto primal_x = mesh.triangles[coparent0y][k];
      auto primal_y = mesh.triangles[coparent0y][(k + 1) % 3];
      has_been_crossed[coparent0y] = true;
      has_been_crossed[coparent1y] = true;
      if (graph[primal_x].is_in_restricted && graph[primal_y].is_in_restricted)
        return vec2i{primal_x, primal_y};
    }

    coparent0x = coparent1x;
    coparent1x = cotree[coparent1x].neighbors[0];
    if (coparent1x == -1)
      stopx = true;

    coparent0y = coparent1y;
    coparent1y = cotree[coparent1y].neighbors[0];
    if (coparent1y == -1)
      stopy = true;

    if (stopy && stopx)
      return vec2i{-1, -1};

    if (has_been_crossed[coparent1x] || has_been_crossed[coparent1y])
      return vec2i{-1, -1};
  }
}

vector<int> create_link(const vector<node> &cut_graph, const vec2i &generators,
                        const vector<int> &cut_points) {

  vector<int> left_side = {generators.x};
  vector<int> right_side = {generators.y};
  while (cut_points[left_side.back()] != 1) {
    auto parent_entry = cut_graph[left_side.back()].neighbors[0];
    left_side.push_back(parent_entry);
  }

  while (cut_points[right_side.back()] != 1) {
    auto parent_entry = cut_graph[right_side.back()].neighbors[0];
    right_side.push_back(parent_entry);
  }

  reverse(left_side.begin(), left_side.end());
  left_side.insert(left_side.end(), right_side.begin(), right_side.end());

  return left_side;
}
vector<int> create_link(const vector<node> &cut_graph,
                        const vec2i &generators) {

  vector<int> left_side = {generators.x};
  vector<int> right_side = {generators.y};
  while (!cut_graph[left_side.back()].is_in_restricted) {
    auto parent_entry = cut_graph[left_side.back()].neighbors[0];
    left_side.push_back(parent_entry);
  }

  while (!cut_graph[right_side.back()].is_in_restricted) {
    auto parent_entry = cut_graph[right_side.back()].neighbors[0];
    right_side.push_back(parent_entry);
  }

  reverse(left_side.begin(), left_side.end());
  left_side.insert(left_side.end(), right_side.begin(), right_side.end());

  return left_side;
}

vector<int> find_generators_in_pruned_tree(const shape_data &mesh,
                                           const shape_topology &topology,
                                           const geodesic_solver &solver,
                                           const vector<node> &cut_graph,
                                           const vector<node> &cut_cotree,
                                           const vec2i &gen) {

  auto entry = node_is_adjacent(solver, gen.x, gen.y);
  auto tid0 = topology.v2t[gen.x][entry / 2];
  auto k = find(mesh.triangles[tid0], gen.x);
  auto tid1 = topology.adjacencies[tid0][k];
  auto gen_in_tree = compute_parent_along_cotree(mesh, topology, cut_cotree,
                                                 cut_graph, vec2i{tid0, tid1});

  if (gen_in_tree.x == -1)
    return create_link(cut_graph, gen);

  return {gen_in_tree.x, gen_in_tree.y};
}
vector<int> find_generators_in_pruned_tree(const shape_data &mesh,
                                           const shape_topology &topology,
                                           const geodesic_solver &solver,
                                           const vector<node> &cut_graph,
                                           const vector<node> &cotree,
                                           const vec2i &gen, const bool &fast) {

  if (fast)
    return find_generators_in_pruned_tree(mesh, topology, solver, cut_graph,
                                          cotree, gen);

  auto entry = node_is_adjacent(solver, gen.x, gen.y);
  auto tid0 = topology.v2t[gen.x][entry / 2];
  auto k = find(mesh.triangles[tid0], gen.x);
  auto tid1 = topology.adjacencies[tid0][k];
  auto parent = compute_parent_along_cotree(mesh, topology, cotree, cut_graph,
                                            {tid0, tid1});
  auto coparent0 = parent[tid1];
  auto coparent1 = parent[coparent0];
  auto gen_in_tree = vec2i{-1, -1};
  bool found = false;
  vector<vec2i> primal_generators;
  while (true) {
    auto k = find(topology.adjacencies[coparent0], coparent1);
    auto primal_x = mesh.triangles[coparent0][k];
    auto primal_y = mesh.triangles[coparent0][(k + 1) % 3];
    primal_generators.push_back(vec2i{primal_x, primal_y});
    if (cut_graph[primal_x].is_in_restricted &&
        cut_graph[primal_y].is_in_restricted) {
      gen_in_tree = {primal_x, primal_y};
      found = true;
      break;
    }
    coparent0 = coparent1;
    coparent1 = parent[coparent1];
    if (coparent1 == -1) {
      break;
    }
  }

  auto lambda = flt_max;
  if (!found) {
    for (auto i = 0; i < primal_generators.size(); ++i) {

      auto parent = primal_generators[i].x;
      auto value = cut_graph[parent].lap;
      while (!cut_graph[parent].is_in_restricted) {
        parent = cut_graph[parent].neighbors[0];
        if (parent == -1)
          break;
        value += cut_graph[parent].lap;
        if (value > lambda)
          break;
      }
      if (value > lambda)
        continue;
      parent = primal_generators[i].y;
      value += cut_graph[parent].lap;
      while (!cut_graph[parent].is_in_restricted) {
        parent = cut_graph[parent].neighbors[0];
        if (parent == -1)
          break;
        value += cut_graph[parent].lap;
        if (value > lambda)
          break;
      }
      if (value < lambda) {
        lambda = value;
        gen_in_tree = primal_generators[i];
      }
    }

    return create_link(cut_graph, gen_in_tree);
  }

  return {gen_in_tree.x, gen_in_tree.y};
}

vector<vector<int>>
find_generators(const shape_data &mesh, const shape_topology &topology,
                const geodesic_solver &solver,
                const vector<node> &cut_graph_mixed, const vector<node> &cotree,
                const vector<vec2i> &generators, const bool fast) {
  auto links = vector<vector<int>>(generators.size());
  for (auto i = 0; i < generators.size(); ++i) {
    if (cut_graph_mixed[generators[i].x].is_in_restricted &&
        cut_graph_mixed[generators[i].y].is_in_restricted)
      links[i] = {generators[i].x, generators[i].y};
    else {

      links[i] = find_generators_in_pruned_tree(
          mesh, topology, solver, cut_graph_mixed, cotree, generators[i], fast);
    }
  }

  return links;
}
bool node_is_an_original_leaf(const vector<node> &cut_graph,
                              const int entry_of_node) {
  if (cut_graph[entry_of_node].neighbors[0] == -1)
    return false;
  if (cut_graph[entry_of_node].neighbors.size() == 1)
    return true;

  return false;
}

vector<node> find_original_leaves(const vector<node> &cut_graph) {
  auto leaves = vector<node>{};
  for (auto i = 0; i < cut_graph.size(); ++i) {
    if (node_is_an_original_leaf(cut_graph, i))
      leaves.push_back(cut_graph[i]);
  }
  return leaves;
}
bool node_is_a_leaf_in_cut_tree(const vector<node> &cut_graph,
                                const int entry_of_parent) {
  if (cut_graph[entry_of_parent].neighbors[0] == -1)
    return false;
  if (cut_graph[entry_of_parent].neighbors.size() == 1)
    return true;

  for (auto i = 0; i < cut_graph[entry_of_parent].neighbors.size(); ++i) {
    if (i == 0)
      continue;
    auto child = cut_graph[entry_of_parent].neighbors[i];
    if (!cut_graph[child].is_cut)
      return false;
  }

  return true;
}
bool node_is_a_leaf(const vector<node> &cut_graph, const int entry_of_parent) {

  if (cut_graph[entry_of_parent].neighbors.size() == 1)
    return true;
  if (cut_graph[entry_of_parent].neighbors[0] == -1)
    return false;
  for (auto i = 1; i < cut_graph[entry_of_parent].neighbors.size(); ++i) {
    if (!cut_graph[cut_graph[entry_of_parent].neighbors[i]].is_cut)
      return false;
  }

  return true;
}
bool node_is_a_leaf_in_restricted(const vector<node> &cut_graph,
                                  const int entry_of_node,
                                  bool count_cut_ones = false) {

  if (cut_graph[entry_of_node].neighbors.size() == 1)
    return true;

  for (auto i = 0; i < cut_graph[entry_of_node].neighbors.size(); ++i) {
    if (i == 0)
      continue;
    auto child = cut_graph[entry_of_node].neighbors[i];
    if (!count_cut_ones) {

      if (cut_graph[child].is_in_restricted && !cut_graph[child].is_cut)
        return false;
    } else {

      if (cut_graph[child].is_in_restricted)
        return false;
    }
  }

  return true;
}

vector<int> find_cut_points(const vector<node> &cut_graph, const int size) {
  vector<int> cut_points(size, 0);
  for (auto node : cut_graph) {
    if (node.is_in_restricted)
      cut_points[node.vid] = 1;
  }

  return cut_points;
}

bool are_relatives(const vector<node> &cut_graph, const int entry0,
                   const int entry1, const bool &strict = false) {

  auto nei = cut_graph[entry0].neighbors;
  if (std::find(nei.begin(), nei.end(), entry1) != nei.end())
    return true;

  if (!strict) {
    if (nei[0] != -1) {
      auto children = cut_graph[nei[0]].neighbors;
      for (auto i = 1; i < children.size(); ++i) {
        if (children[i] == entry1)
          return true;
      }
    }
  }
  return false;
}

void make_branch_unprunable(vector<node> &cut_graph, const int vid) {
  cut_graph[vid].unprunable = true;
  cut_graph[vid].is_in_restricted = true;
  cut_graph[vid].is_cut = false;
  auto parent_entry = cut_graph[vid].neighbors[0];
  while (parent_entry != -1) {
    cut_graph[parent_entry].unprunable = true;
    cut_graph[parent_entry].is_in_restricted = true;
    cut_graph[parent_entry].is_cut = false;
    parent_entry = cut_graph[parent_entry].neighbors[0];
  }
}

void update_tree(const int vid0, const int vid1, vector<node> &cut_graph) {

  make_branch_unprunable(cut_graph, vid0);
  make_branch_unprunable(cut_graph, vid1);
}

void update_tree(const vector<vector<int>> &links,
                 const vector<pair<float, int>> &lap, vector<node> &cut_graph) {
  auto count_added = 0;
  auto added = vector<int>{};
  for (auto i = 0; i < links.size(); ++i) {
    if (links[i][0] == -1)
      continue;

    ++count_added;
    for (auto j = 0; j < links[i].size(); ++j) {

      auto entry = links[i][j];
      auto prev_entry = (j > 0) ? links[i][j - 1] : -1;
      auto next_entry = (j < links[i].size() - 1) ? links[i][j + 1] : -1;
      cut_graph[entry].unprunable = true;
      cut_graph[entry].is_in_restricted = true;
      cut_graph[entry].is_cut = false;

      if (prev_entry >= 0 && !are_relatives(cut_graph, entry, prev_entry, true))
        cut_graph[entry].neighbors.push_back(prev_entry);

      if (next_entry >= 0 && !are_relatives(cut_graph, entry, next_entry, true))
        cut_graph[entry].neighbors.push_back(next_entry);
    }

    update_tree(links[i][0], links[i].back(), cut_graph);
  }
  if (links.size() > 0) {
    std::cout << "" << std::endl;
    printf("%d generator found out of %d", count_added, (int)links.size());
    std::cout << "" << std::endl;
  }
}

vector<node> find_leaves(const vector<node> &cut_graph) {
  vector<node> leaves;
  for (auto i = 0; i < cut_graph.size(); ++i) {
    if (cut_graph[i].is_cut)
      continue;
    if (node_is_a_leaf(cut_graph, i))
      leaves.push_back(cut_graph[i]);
  }
  return leaves;
}

vector<node> find_leaves_in_restricted(const vector<node> &cut_graph) {
  auto leaves = vector<node>{};
  for (auto i = 0; i < cut_graph.size(); ++i) {
    if (cut_graph[i].is_cut)
      continue;
    if (!cut_graph[i].is_in_restricted)
      continue;
    if (node_is_a_leaf_in_restricted(cut_graph, i))
      leaves.push_back(cut_graph[i]);
  }
  return leaves;
}

int count_children(const vector<node> &cut_graph, const int parent_entry) {
  auto count = 0;

  for (auto nei : cut_graph[parent_entry].neighbors) {
    if (nei != -1 && !cut_graph[nei].is_cut)
      ++count;
  }
  return count;
}

int count_children_in_restricted(const vector<node> &cut_graph,
                                 const int parent_entry) {

  auto count = 0;
  for (auto nei : cut_graph[parent_entry].neighbors) {
    if (nei != -1 && cut_graph[nei].is_in_restricted && !cut_graph[nei].is_cut)
      ++count;
  }

  return count;
}

bool is_a_branching_node(const vector<node> &cut_graph, const int entry) {
  auto valence = count_children(cut_graph, entry);
  return (valence >= 3);
}

bool is_a_branching_node_in_restricted(const vector<node> &cut_graph,
                                       const int entry) {
  auto valence = count_children_in_restricted(cut_graph, entry);
  return (valence >= 3);
}

int find_parent(const vector<bool> &pushed,
                vector<geodesic_solver::graph_edge> &nbr,
                const vector<pair<float, int>> &hess) {
  auto candidate = vector<pair<float, int>>{};

  for (auto i = 0; i < nbr.size(); ++i) {
    if (i % 2 == 1)
      continue;
    if (pushed[nbr[i].node]) {
      candidate.push_back(hess[nbr[i].node]);
    }
  }

  if (candidate.size() == 0)
    return -1;
  else
    sort(candidate.begin(), candidate.end());

  return candidate[0].second;
}

int find_parent(const vector<bool> &pushed,
                const vector<geodesic_solver::graph_edge> &nbr,
                const vector<int> &cut_points,
                const vector<pair<float, int>> &lap,
                const vector<float> &distances) {
  auto red_candidate = vector<pair<float, int>>{};
  auto white_candidate = vector<pair<float, int>>{};
  for (auto i = 0; i < nbr.size(); ++i) {
    if (i % 2 == 1)
      continue;
    if (pushed[nbr[i].node]) {
      if (cut_points[nbr[i].node] == 1)
        red_candidate.push_back(lap[nbr[i].node]);
      else
        white_candidate.push_back({distances[nbr[i].node], nbr[i].node});
    }
  }

  if (red_candidate.size() == 0 && white_candidate.size() == 0)
    return -1;
  else if (red_candidate.size() > 0) {
    sort(red_candidate.begin(), red_candidate.end());
    return red_candidate[0].second;
  } else {
    sort(white_candidate.begin(), white_candidate.end());
    return white_candidate.back().second;
  }
}

bool is_in_restricted(const node &curr_leave, const float &grad_threshold,
                      const float &lap_threshold, const bool &use_grad,
                      const bool &use_lap) {
  if (!use_lap && !use_grad)
    return true;

  if (use_grad && use_lap) {
    if (curr_leave.angle_deviation < grad_threshold ||
        curr_leave.lap > lap_threshold)
      return false;
  } else if (use_grad) {
    if (curr_leave.angle_deviation < grad_threshold)

      return false;

  } else {
    if (curr_leave.lap > lap_threshold)
      return false;
  }

  return true;
}
// optimize
void update_restricted(const float &grad_threshold, const float &lap_threshold,
                       const bool &use_grad, const bool &use_lap,
                       const bool &propagate, const int entry_of_min,
                       vector<node> &cut_graph) {
  for (auto i = 0; i < cut_graph.size(); ++i) {
    if (cut_graph[i].unprunable)
      continue;

    cut_graph[i].is_in_restricted = false;
  }
  auto Q = std::deque<node>{};
  auto has_been_pushed = vector<bool>(cut_graph.size(), false);
  if (propagate)
    Q = {cut_graph[entry_of_min]};
  else {
    auto leaves = find_original_leaves(cut_graph);
    for (auto leaf : leaves) {
      Q.push_back(leaf);
    }
  }

  while (!Q.empty()) {
    auto leaf = Q.back();
    Q.pop_back();
    if (has_been_pushed[leaf.vid])
      continue;

    has_been_pushed[leaf.vid] = true;
    if (is_in_restricted(leaf, grad_threshold, lap_threshold, use_grad,
                         use_lap)) {
      cut_graph[leaf.vid].is_in_restricted = true;

      if (propagate) {
        for (auto i = 0; i < cut_graph[leaf.vid].neighbors.size(); ++i) {
          if (i == 0)
            continue;
          if (!cut_graph[cut_graph[leaf.vid].neighbors[i]].is_cut)
            Q.push_back(cut_graph[cut_graph[leaf.vid].neighbors[i]]);
        }
      } else {
        auto parent = leaf.neighbors[0];
        while (parent != -1) {
          cut_graph[parent].is_in_restricted = true;
          cut_graph[parent].is_cut = false;
          has_been_pushed[cut_graph[parent].vid] = true;
          parent = cut_graph[parent].neighbors[0];
        }
      }
    } else {
      if (!propagate) {
        auto parent = leaf.neighbors[0];
        if (parent != -1) {
          if (node_is_a_leaf_in_restricted(cut_graph, parent))
            Q.push_back(cut_graph[parent]);
        }
      }
    }
  }
}
bool update_node(const shape_data &mesh, const shape_topology &topology,
                 const geodesic_solver &solver,
                 const vector<pair<float, int>> &lap, const int child,
                 const int parent, vector<node> &cut_graph,
                 vector<vec3i> &edge_mask) {
  if (cut_graph[child].vid == -1) {
    auto curr_node = node{};
    curr_node.vid = child;
    curr_node.lap = lap[child].first;
    curr_node.neighbors = {parent};
    curr_node.is_in_restricted = cut_graph[child].is_in_restricted;
    auto entry = node_is_adjacent(solver, parent, child);
    auto tid0 = topology.v2t[parent][entry / 2];
    auto k = find(mesh.triangles[tid0], parent);
    auto tid1 = opposite_face(mesh.triangles, topology.adjacencies, tid0,
                              mesh.triangles[tid0][(k + 2) % 3]);
    auto h0 = find(topology.adjacencies[tid0], tid1);
    auto h1 = find(topology.adjacencies[tid1], tid0);
    edge_mask[tid0][h0] = 1;
    edge_mask[tid1][h1] = 1;
    cut_graph[child] = curr_node;
    return true;
  }
  return false;
}
void mixed_spanning_tree(const shape_data &mesh, const shape_topology &topology,
                         const geodesic_solver &solver,
                         const vector<node> &lap_graph,
                         const vector<float> &distances,
                         const vector<pair<float, int>> &lap,
                         const float &max_lap, const int entry_of_min,
                         vector<node> &mixed_graph, vector<vec3i> &edge_mask) {
  time_function();
  auto priority_for_cut_locus =
      [](const float &max_lap, const float &max_dist, const float &lap_at_vid)

  { return -lap_at_vid + max_dist + max_lap; };

  typedef pair<float, int> hess_at_v;
  std::priority_queue<hess_at_v, vector<hess_at_v>, std::less<hess_at_v>> Q;

  mixed_graph.resize(mesh.positions.size(), node{});
  edge_mask.resize(mesh.triangles.size(), zero3i);
  auto max_dist = distances[entry_of_min];
  auto temp_node = node{};
  temp_node.vid = entry_of_min;
  temp_node.neighbors = {-1};
  temp_node.lap = lap[entry_of_min].first;
  mixed_graph[entry_of_min] = temp_node;
  for (auto i = 0; i < mesh.positions.size(); ++i) {
    if (lap_graph[i].is_in_restricted)
      mixed_graph[i].is_in_restricted = true;
  }
  Q.push(std::make_pair(
      priority_for_cut_locus(max_lap, max_dist, lap[entry_of_min].first),
      entry_of_min));
  while (!Q.empty()) {
    auto root = Q.top();
    Q.pop();
    auto nbr = solver.graph[root.second];
    for (auto i = 0; i < solver.graph[root.second].size(); ++i) {
      auto nei = solver.graph[root.second][i];
      if (i % 2 == 1 || !update_node(mesh, topology, solver, lap, nei.node,
                                     root.second, mixed_graph, edge_mask))
        continue;
      mixed_graph[root.second].neighbors.push_back(nei.node);
      if (mixed_graph[nei.node].is_in_restricted) {
        Q.push({priority_for_cut_locus(max_lap, max_dist, lap[nei.node].first),
                nei.node});
      } else {
        Q.push({distances[nei.node], nei.node});
      }
    }
  }
}

inline bool is_in_tree(const vector<node> &graph, int vid) {
  return (graph[vid].vid != -1);
}
void spanning_tree(const shape_data &mesh, const shape_topology &topology,
                   const geodesic_solver &solver,
                   const vector<pair<float, int>> &lap,
                   const vector<vec3f> &gradients, const int entry_of_root,
                   vector<node> &graph) {
  typedef pair<float, int> hess_at_v;
  std::priority_queue<hess_at_v, vector<hess_at_v>, std::greater<hess_at_v>> Q;
  auto temp_node = node{};
  temp_node.vid = entry_of_root;
  temp_node.lap = lap[entry_of_root].first;
  temp_node.neighbors = {-1};
  temp_node.angle_deviation = graph[entry_of_root].angle_deviation;
  auto nbr = solver.graph[temp_node.vid];
  graph[entry_of_root] = temp_node;
  Q.push(lap[entry_of_root]);

  while (!Q.empty()) {
    auto curr = Q.top();
    Q.pop();
    auto parent = curr.second;
    auto nbr = solver.graph[parent];
    for (auto i = 0; i < nbr.size(); ++i) {
      if (i % 2 == 1 || is_in_tree(graph, nbr[i].node))
        continue;
      // add child to tree
      auto child = nbr[i].node;
      auto temp_node = node{};
      temp_node.vid = child;
      temp_node.lap = lap[child].first;
      temp_node.angle_deviation = graph[child].angle_deviation;
      temp_node.neighbors = {parent};
      graph[child] = temp_node;
      graph[parent].neighbors.push_back(child);
      Q.push(lap[child]);
    }
  }
}

float field_on_tid(const shape_data &mesh, const vector<float> &field,
                   const int tid) {

  return 0.33 * (field[mesh.triangles[tid].x] + field[mesh.triangles[tid].y] +
                 field[mesh.triangles[tid].z]);
}

int find_coparent(const shape_topology &topology, const vector<node> &tree,
                  const vector<bool> &has_been_pushed,
                  const vector<int> &mapping,
                  const vector<vec3i> &crossed_edges, const int &tid,
                  const bool greater) {
  auto lambda = (greater) ? flt_min : flt_max;
  auto candidate = -1;
  if (greater) {
    for (auto i = 0; i < 3; ++i) {
      auto nei = topology.adjacencies[tid][i];
      if (has_been_pushed[nei] && crossed_edges[tid][i] == 0) {
        auto value = tree[mapping.at(nei)].lap;
        if (value > lambda) {
          lambda = value;
          candidate = nei;
        }
      }
    }
  } else {
    for (auto i = 0; i < 3; ++i) {
      auto nei = topology.adjacencies[tid][i];
      if (has_been_pushed[nei] && crossed_edges[tid][i] == 0) {
        auto value = tree[mapping[nei]].lap;
        if (value < lambda) {
          lambda = value;
          candidate = nei;
        }
      }
    }
  }

  return candidate;
}

bool update_conode(const shape_data &mesh, const shape_topology &topology,
                   const vector<float> &distances, const int child,
                   const int parent, vector<node> &cut_graph,
                   vector<vec3i> &edge_mask) {
  if (cut_graph[child].vid == -1) {
    auto curr_node = node{};
    curr_node.vid = child;
    curr_node.lap = field_on_tid(mesh, distances, child);
    curr_node.neighbors = {parent};
    auto k = find(topology.adjacencies[child], parent);
    edge_mask[child][k] = -1;
    k = find(topology.adjacencies[parent], child);
    edge_mask[parent][k] = -1;
    cut_graph[child] = curr_node;
    return true;
  }
  return false;
}

void spanning_cotree(const shape_data &mesh, const shape_topology &topology,
                     const vector<float> &distances, const mesh_point &source,
                     vector<node> &graph, vector<vec3i> &edge_mask) {
  time_function();
  auto seed = source.face;
  auto curr_node = node{};
  graph.resize(mesh.triangles.size(), curr_node);
  curr_node.vid = seed;
  curr_node.lap = field_on_tid(mesh, distances, seed);
  curr_node.neighbors = {-1};
  typedef pair<float, int> field_at_v;
  std::priority_queue<field_at_v, vector<field_at_v>, std::greater<field_at_v>>
      Q;
  Q.push(std::make_pair(field_on_tid(mesh, distances, seed), seed));
  while (!Q.empty()) {
    auto curr = Q.top();
    Q.pop();

    for (auto i = 0; i < 3; ++i) {
      auto nei = topology.adjacencies[curr.second][i];
      if (edge_mask[curr.second][i] == 0 &&
          update_conode(mesh, topology, distances, nei, curr.second, graph,
                        edge_mask)) {
        graph[curr.second].neighbors.push_back(nei);
        Q.push(std::make_pair(field_on_tid(mesh, distances, nei), nei));
      }
    }
  }
}

bool prune_branch_from_top(vector<node> &cut_graph, const int vid) {
  std::deque<int> Q;
  if (cut_graph[vid].unprunable)
    return false;

  auto branch = vector<int>{vid};
  for (auto i = 0; i < cut_graph[vid].neighbors.size(); ++i) {
    if (i == 0)
      continue;

    Q.push_back(cut_graph[vid].neighbors[i]);
  }

  while (!Q.empty()) {
    auto curr = Q.back();
    Q.pop_back();
    if (cut_graph[curr].unprunable)
      return false;
    branch.push_back(curr);
    for (auto j = 0; j < cut_graph[curr].neighbors.size(); ++j) {
      if (j == 0)
        continue;

      Q.push_back(cut_graph[curr].neighbors[j]);
    }
  }

  for (auto child : branch) {
    cut_graph[child].is_in_restricted = false;
  }

  return true;
}
bool prune_branch_from_top_in_restricted(vector<node> &cut_graph,
                                         const int vid) {
  std::deque<int> Q;
  if (cut_graph[vid].unprunable)
    return false;

  auto branch = vector<int>{vid};

  for (auto i = 0; i < cut_graph[vid].neighbors.size(); ++i) {
    if (i == 0)
      continue;

    Q.push_back(cut_graph[vid].neighbors[i]);
  }

  while (!Q.empty()) {
    auto curr = Q.back();
    Q.pop_back();
    if (cut_graph[curr].unprunable)
      return false;
    branch.push_back(curr);
    for (auto j = 0; j < cut_graph[curr].neighbors.size(); ++j) {
      if (j == 0)
        continue;
      if (cut_graph[cut_graph[curr].neighbors[j]].is_in_restricted &&
          !cut_graph[cut_graph[curr].neighbors[j]].is_cut)
        Q.push_back(cut_graph[curr].neighbors[j]);
    }
  }

  for (auto child : branch) {
    cut_graph[child].is_in_restricted = false;
  }

  return true;
}

void cut_short_branches(vector<node> &cut_graph, const int depth) {
  auto leaves = find_leaves(cut_graph);

  for (auto leaf : leaves) {
    auto curr = leaf;
    auto count_depth = 0;

    while (count_depth < depth) {
      auto parent = curr.neighbors[0];
      if (parent == -1) {

        break;
      }
      if (is_a_branching_node(cut_graph, parent)) {
        prune_branch_from_top(cut_graph, curr.vid);
        break;
      }
      ++count_depth;

      curr = cut_graph[parent];
    }
  }
}

void cut_short_branches_in_restricted(vector<node> &cut_graph,
                                      const int depth) {
  // update_depth(curr_tree);
  auto leaves = find_leaves_in_restricted(cut_graph);
  for (auto leaf : leaves) {
    auto curr = leaf;
    auto count_depth = 0;

    while (count_depth < depth) {
      auto parent = curr.neighbors[0];
      if (parent == -1) {

        break;
      }
      if (is_a_branching_node_in_restricted(cut_graph, parent)) {
        prune_branch_from_top_in_restricted(cut_graph, curr.vid);
        break;
      }
      ++count_depth;

      curr = cut_graph[parent];
    }
  }
}

vector<int> find_relatives(const vector<node> &cut_graph, const int entry) {
  auto relatives = vector<int>{};

  auto neis = cut_graph[entry].neighbors;
  for (auto nei : neis) {
    if (nei == -1)
      continue;
    if (!cut_graph[nei].is_cut && cut_graph[nei].is_in_restricted)
      relatives.push_back(nei);
  }

  return relatives;
}

void uncut_tree(vector<node> &cut_graph) {
  for (auto i = 0; i < cut_graph.size(); ++i) {
    cut_graph[i].is_cut = false;
  }
}

void uncut_tree_in_restricted(vector<node> &cut_graph) {
  for (auto i = 0; i < cut_graph.size(); ++i) {
    if (!cut_graph[i].is_in_restricted)
      continue;
    cut_graph[i].is_cut = false;
  }
}

void thinning(const shape_data &mesh, const geodesic_solver &solver,
              const vector<pair<float, int>> &lap, vector<node> &cut_graph) {
  time_function();
  auto leaves = find_original_leaves(cut_graph);
  std::deque<node> Q;
  auto is_in_queue = vector<bool>(mesh.positions.size(), false);
  for (auto leaf : leaves) {

    Q.push_back(leaf);
    is_in_queue[leaf.vid] = true;
  }

  while (!Q.empty()) {
    auto curr = Q.front();
    Q.pop_front();

    auto nbr = solver.graph[curr.vid];
    for (auto i = 0; i < nbr.size(); ++i) {
      if (i % 2)
        continue;
      if (cut_graph[nbr[i].node].is_cut)
        continue;
      if (!are_relatives(cut_graph, curr.vid, nbr[i].node)) {

        if (node_is_a_leaf_in_cut_tree(cut_graph, nbr[i].node)) {
          if (lap[curr.vid].first > lap[nbr[i].node].first &&
              !cut_graph[curr.vid].unprunable) {
            cut_graph[curr.vid].is_cut = true;
            cut_graph[curr.vid].is_in_restricted = false;
            auto parent = cut_graph[curr.vid].neighbors[0];
            if (parent != -1) {
              if (node_is_a_leaf_in_cut_tree(cut_graph, parent) &&
                  !is_in_queue[parent] && !cut_graph[parent].unprunable) {
                Q.push_back(cut_graph[parent]);
                is_in_queue[parent] = true;
              }
            }
          }
        } else {
          cut_graph[curr.vid].is_cut = true;
          cut_graph[curr.vid].is_in_restricted = false;
          auto parent = cut_graph[curr.vid].neighbors[0];
          if (parent != -1) {
            if (node_is_a_leaf_in_cut_tree(cut_graph, parent) &&
                !is_in_queue[parent] && !cut_graph[parent].unprunable) {
              Q.push_back(cut_graph[parent]);
              is_in_queue[parent] = true;
            }
          }
        }
      }
    }
    is_in_queue[curr.vid] = false;
  }
}

void thinning_in_restricted(const shape_data &mesh,
                            const geodesic_solver &solver,
                            const vector<pair<float, int>> &lap,
                            vector<node> &cut_graph) {
  uncut_tree_in_restricted(cut_graph);
  auto leaves = find_leaves_in_restricted(cut_graph);
  std::deque<node> Q;
  auto is_in_queue = vector<bool>(mesh.positions.size(), false);
  for (auto leaf : leaves) {

    Q.push_back(leaf);
    is_in_queue[leaf.vid] = true;
  }

  while (!Q.empty()) {
    auto curr = Q.front();
    Q.pop_front();

    auto nbr = solver.graph[curr.vid];
    for (auto i = 0; i < nbr.size(); ++i) {
      if (i % 2)
        continue;

      if (!are_relatives(cut_graph, curr.vid, nbr[i].node)) {

        if (node_is_a_leaf_in_restricted(cut_graph, nbr[i].node)) {
          if (lap[curr.vid].first > lap[nbr[i].node].first &&
              !cut_graph[nbr[i].node].unprunable) {
            cut_graph[curr.vid].is_cut = true;
            cut_graph[curr.vid].is_in_restricted = false;
            auto parent = cut_graph[curr.vid].neighbors[0];
            if (parent != -1) {
              if (node_is_a_leaf_in_restricted(cut_graph, parent) &&
                  !is_in_queue[parent] && !cut_graph[parent].unprunable) {
                Q.push_back(cut_graph[parent]);
                is_in_queue[parent] = true;
              }
            }
          }
        } else {
          cut_graph[curr.vid].is_cut = true;
          cut_graph[curr.vid].is_in_restricted = false;
          auto parent = cut_graph[curr.vid].neighbors[0];
          if (parent != -1) {
            if (node_is_a_leaf_in_restricted(cut_graph, parent) &&
                !is_in_queue[parent] && !cut_graph[parent].unprunable) {
              Q.push_back(cut_graph[parent]);
              is_in_queue[parent] = true;
            }
          }
        }
      }
    }
    is_in_queue[curr.vid] = false;
  }
}
float connect_centroids(const shape_data &mesh, const shape_topology &topology,
                        const int tid0, const int tid1) {
  auto flat_tid = init_flat_triangle(mesh.positions, mesh.triangles[tid0]);
  auto flat_opp =
      unfold_face(mesh.triangles, mesh.positions, flat_tid, tid0, tid1);

  auto k = find(topology.adjacencies[tid0], tid1);
  auto left = flat_tid[k];
  auto right = flat_tid[(k + 1) % 3];
  auto c0 = interpolate_triangle(flat_tid[0], flat_tid[1], flat_tid[2],
                                 vec2f{0.33, 0.33});
  auto c1 = interpolate_triangle(flat_opp[0], flat_opp[1], flat_opp[2],
                                 vec2f{0.33, 0.33});

  // intersection
  auto a = c1 - c0;
  auto b = right - left;
  auto d = right - c0;
  auto det = a.x * b.y - a.y * b.x;
  assert(det);
  auto lerp = (a.x * d.y - a.y * d.x) / det;
  return lerp;
}

std::tuple<vector<vec2i>, vector<vec3f>, vector<vec2i>, vector<vec3f>>
compute_splitted_branches(const shape_data &mesh,
                          const vector<node> &cut_graph) {

  auto branches = vector<vec2i>{};
  auto node_pos = vector<vec3f>{};
  auto restricted_branches = vector<vec2i>{};
  auto restricted_node_pos = vector<vec3f>{};
  auto map = unordered_map<int, int>{};
  auto restr_map = unordered_map<int, int>{};
  for (auto i = 0; i < cut_graph.size(); ++i) {
    if (cut_graph[i].is_cut)
      continue;
    auto n = mesh.normals[cut_graph[i].vid];
    if (cut_graph[i].is_in_restricted) {

      restricted_node_pos.push_back(mesh.positions[i] + 0.0001 * n);
      restr_map[i] = (int)(restricted_node_pos.size() - 1);
    } else {

      node_pos.push_back(mesh.positions[i] + 0.0001 * n);
      map[i] = (int)(node_pos.size() - 1);
    }
  }
  for (auto j = 0; j < cut_graph.size(); ++j) {
    if (cut_graph[j].is_cut)
      continue;

    for (auto i = 0; i < cut_graph[j].neighbors.size(); ++i) {
      if (i == 0)
        continue;
      auto child_entry = cut_graph[j].neighbors[i];
      if (cut_graph[child_entry].is_cut)
        continue;
      if (cut_graph[j].is_in_restricted &&
          cut_graph[child_entry].is_in_restricted) {
        restricted_branches.push_back(
            {restr_map.at(j), restr_map.at(child_entry)});
      } else if (cut_graph[j].is_in_restricted) {
        auto size = (int)(restricted_node_pos.size());
        restricted_node_pos.resize(size + 1);
        restricted_node_pos[size] = node_pos[map.at(child_entry)];
        restricted_branches.push_back({restr_map.at(j), size});
      } else if (cut_graph[child_entry].is_in_restricted) {
        auto size = (int)(node_pos.size());
        node_pos.resize(size + 1);
        node_pos[size] = restricted_node_pos[restr_map.at(child_entry)];
        branches.push_back({map.at(j), size});
      } else
        branches.push_back({map.at(j), map.at(child_entry)});
    }
  }
  return {branches, node_pos, restricted_branches, restricted_node_pos};
}

std::pair<vector<vec2i>, vector<vec3f>>
compute_branches(const shape_data &mesh, const vector<node> &cut_graph) {

  auto branches = vector<vec2i>{};
  auto node_pos = vector<vec3f>(cut_graph.size());
  for (auto i = 0; i < cut_graph.size(); ++i) {
    if (cut_graph[i].is_cut)
      continue;
    if (!cut_graph[i].is_in_restricted)
      continue;
    auto n = mesh.normals[cut_graph[i].vid];
    node_pos[i] = mesh.positions[cut_graph[i].vid] + 0.0001 * n;
    for (auto j = 0; j < cut_graph[i].neighbors.size(); ++j) {
      if (j == 0)
        continue;
      auto child_entry = cut_graph[i].neighbors[j];
      if (cut_graph[child_entry].is_cut)
        continue;
      if (!cut_graph[child_entry].is_in_restricted)
        continue;
      branches.push_back({i, child_entry});
    }
  }
  return {branches, node_pos};
}

std::pair<vector<vec2i>, vector<vec3f>>
compute_smoothed_branches(const shape_data &mesh,
                          const shape_topology &topology,
                          const vector<node> &cut_graph,
                          const vector<mesh_point> &new_one, const bool &draw) {
  auto branches = vector<vec2i>{};
  auto node_pos = vector<vec3f>{};
  auto map = std::unordered_map<int, int>{};
  for (auto i = 0; i < cut_graph.size(); ++i) {
    if (cut_graph[i].is_cut || !cut_graph[i].is_in_restricted)
      continue;
    auto n =
        (draw) ? tid_normal(mesh.triangles, mesh.positions, new_one.at(i).face)
               : zero3f;
    node_pos.push_back(eval_position(mesh, new_one.at(i)) + 0.001 * n);

    map[i] = (int)(node_pos.size() - 1);
  }
  for (auto i = 0; i < cut_graph.size(); ++i) {
    if (cut_graph[i].is_cut || !cut_graph[i].is_in_restricted)
      continue;
    for (auto h = 1; h < cut_graph[i].neighbors.size(); ++h) {

      auto child_entry = cut_graph[i].neighbors[h];
      if (cut_graph[child_entry].is_cut ||
          !cut_graph[child_entry].is_in_restricted)
        continue;

      branches.push_back({map.at(i), map.at(child_entry)});
    }
  }
  return {branches, node_pos};
}

std::pair<vector<vec2i>, vector<vec3f>>
compute_cobranches(const shape_data &mesh, const shape_topology &topology,
                   const vector<node> &cut_graph) {

  auto branches = vector<vec2i>{};
  auto node_pos = vector<vec3f>(cut_graph.size());
  for (auto i = 0; i < cut_graph.size(); ++i) {
    auto t_n = tid_normal(mesh.triangles, mesh.positions, cut_graph[i].vid);
    node_pos[i] =
        tid_centroid(mesh.triangles, mesh.positions, cut_graph[i].vid) +
        0.0001 * t_n;
    for (auto j = 0; j < cut_graph[i].neighbors.size(); ++j) {
      if (j == 0)
        continue;
      auto child_entry = cut_graph[i].neighbors[j];
      auto size = (int)node_pos.size();
      node_pos.resize(size + 1);

      auto lerp = connect_centroids(mesh, topology, cut_graph[i].vid,
                                    cut_graph[child_entry].vid);

      auto k = find(topology.adjacencies[cut_graph[i].vid],
                    cut_graph[child_entry].vid);
      auto e_n =
          lerp * mesh.normals[mesh.triangles[cut_graph[i].vid][k]] +
          (1 - lerp) *
              mesh.normals[mesh.triangles[cut_graph[i].vid][(k + 1) % 3]];
      node_pos[size] =
          lerp * mesh.positions[mesh.triangles[cut_graph[i].vid][k]] +
          (1 - lerp) *
              mesh.positions[mesh.triangles[cut_graph[i].vid][(k + 1) % 3]] +
          0.0001 * e_n;
      branches.push_back({i, size});
      branches.push_back({size, child_entry});
    }
  }
  return {branches, node_pos};
}

vector<vec2i> export_cobranches(const shape_data &mesh,
                                const vector<node> &cut_graph) {
  auto branches = vector<vec2i>{};
  for (auto i = 0; i < cut_graph.size(); ++i) {
    auto tid = cut_graph[i].vid;
    for (auto j = 0; j < cut_graph[i].neighbors.size(); ++j) {
      if (j == 0)
        continue;
      auto child_entry = cut_graph[i].neighbors[j];
      auto eid = common_edge(mesh.triangles[tid],
                             mesh.triangles[cut_graph[child_entry].vid]);
      branches.push_back(eid);
    }
  }
  return branches;
}

void p_transp_along_path(const shape_data &mesh, const shape_topology &topology,
                         const geodesic_solver &solver,
                         const geodesic_path &path, vec3f v) {
  for (auto i = 0; i < path.strip.size() - 1; ++i) {
    parallel_transport(mesh, topology, solver, v, path.strip[i],
                       path.strip[i + 1], T2T);
  }
}

vec3f interpolate_grad(const shape_data &mesh, const shape_topology &topology,
                       const geodesic_solver &solver,
                       const vector<vec3f> &gradients, const mesh_point &p) {
  auto bary = make_bary(p.uv);
  auto gx = gradients[mesh.triangles[p.face].x];
  auto gy = gradients[mesh.triangles[p.face].y];
  auto gz = gradients[mesh.triangles[p.face].z];
  parallel_transport(mesh, topology, solver, gx, mesh.triangles[p.face].x,
                     p.face, V2T);
  parallel_transport(mesh, topology, solver, gy, mesh.triangles[p.face].y,
                     p.face, V2T);
  parallel_transport(mesh, topology, solver, gz, mesh.triangles[p.face].z,
                     p.face, V2T);
  return bary.x * gx + bary.y * gy + bary.z * gz;
}

vec3f dir_of_smoothing(const shape_data &mesh, const shape_topology &topology,
                       const geodesic_solver &solver,
                       const vector<vec3f> &gradients, const mesh_point &prev,
                       const mesh_point &curr, const mesh_point &next,
                       const bool use_gradient) {

  auto prev_pos = eval_position(mesh, prev);
  auto curr_pos = eval_position(mesh, curr);
  auto next_pos = eval_position(mesh, next);
  auto [is_vert, k] = point_is_vert(curr);
  auto n = (is_vert)
               ? mesh.normals[forced_vert_from_point(mesh.triangles, curr)]
               : tid_normal(mesh.triangles, mesh.positions, curr.face);

  auto e1 = normalize(curr_pos - prev_pos);
  auto e2 = normalize(next_pos - curr_pos);
  auto l1 = length(e1);
  auto l2 = length(e2);
  auto normal_dir = project_vec(2 / (l1 + l2) * (e2 - e1), n);

  if (use_gradient) {
    auto grad = interpolate_grad(mesh, topology, solver, gradients, curr);
    return normalize(normal_dir) * dot(normal_dir, grad);
  } else
    return normal_dir;
}
vec3f dir_of_smoothing(const shape_data &mesh, const shape_topology &topology,
                       const geodesic_solver &solver,
                       const dual_geodesic_solver &dual_solver,
                       const vector<vec3f> &gradients, mesh_point &p0,
                       mesh_point &p1, const bool use_gradient,
                       vec3f &reference) {
  auto [p0_is_vert, k0] = point_is_vert(p0);
  auto [p1_is_vert, k1] = point_is_vert(p1);
  auto grad = zero3f;
  if (p0_is_vert && p1_is_vert) {
    auto p0_vid = mesh.triangles[p0.face][k0];
    auto p1_vid = mesh.triangles[p1.face][k1];
    parallel_transport(mesh, topology, solver, reference, p1_vid, p0_vid, V2V);
    reference *= -1;
    if (use_gradient)
      grad = gradients[p0_vid];

  } else {
    auto path =
        compute_geodesic_path(mesh, topology, dual_solver, solver, p1, p0);
    p_transp_along_path(mesh, topology, solver, path, reference);
    reference *= -1;
    if (use_gradient)
      grad = interpolate_grad(mesh, topology, solver, gradients, p0);
  }
  auto dir =
      (use_gradient) ? normalize(reference) * dot(reference, grad) : reference;
  return dir;
}
vec3f average_grad_3D(const shape_data &mesh, const shape_topology &topology,
                      const dual_geodesic_solver &dual_solver,
                      const geodesic_solver &solver,
                      const vector<vec3f> &gradients,
                      const vector<int> &entries,
                      const vector<mesh_point> &mapping, const mesh_point &p0) {

  auto grad = zero3f;

  auto n = tid_normal(mesh.triangles, mesh.positions, p0.face);
  for (auto i = 0; i < entries.size(); ++i) {

    auto curr = mapping[entries[i]];
    auto path =
        compute_geodesic_path(mesh, topology, dual_solver, solver, p0, curr);
    auto len =
        path_length(path, mesh.triangles, mesh.positions, topology.adjacencies);
    auto pos = path_positions(mesh, topology, path);
    auto e1 = pos[1] - pos[0];
    e1 = normalize(project_vec(e1, n));
    e1 *= len;

    grad += e1;
  }

  return grad;
}
vector<mesh_point> init_mapping(const shape_data &mesh,
                                const shape_topology &topology,
                                const vector<node> &cut_graph) {
  vector<mesh_point> mapping(cut_graph.size());
  for (auto i = 0; i < cut_graph.size(); ++i) {
    if (cut_graph[i].is_cut || !cut_graph[i].is_in_restricted)
      mapping[i] = {0, {1, 0}};
    else
      mapping[i] =
          make_point_from_vert(mesh.triangles, topology.v2t, cut_graph[i].vid);
  }

  return mapping;
}
void handle_collapsing_nodes(vector<node> &cut_graph, vector<int> &relatives,
                             const shape_data &mesh,
                             const vector<mesh_point> &mapping,
                             const int curr) {
  auto prev_pos = eval_position(mesh, mapping[relatives[0]]);
  auto curr_pos = eval_position(mesh, mapping[curr]);
  auto next_pos = eval_position(mesh, mapping[relatives[1]]);

  if (length(prev_pos - curr_pos) < 1e-5) {
    auto prev_relatives = find_relatives(cut_graph, relatives[0]);
    if (prev_relatives.size() == 2) {
      auto new_neighbor =
          (prev_relatives[0] == curr) ? prev_relatives[1] : prev_relatives[0];
      cut_graph[relatives[0]].is_cut = true;
      cut_graph[new_neighbor].neighbors.push_back(curr);
      relatives = {new_neighbor, relatives[1]};
    }
  }

  if (length(next_pos - curr_pos) < 1e-5) {
    auto next_relatives = find_relatives(cut_graph, relatives[1]);
    if (next_relatives.size() == 2) {
      auto new_neighbor =
          (next_relatives[0] == curr) ? next_relatives[1] : next_relatives[0];
      cut_graph[relatives[1]].is_cut = true;
      cut_graph[new_neighbor].neighbors.push_back(curr);
      relatives = {relatives[0], new_neighbor};
    }
  }
}

vector<mesh_point>
smoothed_cut_locus(const shape_data &mesh, const shape_topology &topology,
                   const geodesic_solver &solver,
                   const dual_geodesic_solver &dual_solver,
                   vector<node> &cut_graph, const vector<vec3f> &gradients,
                   const int it, float &factor, const bool &use_gradient) {
  auto mapping = init_mapping(mesh, topology, cut_graph);
  auto subsitute = vector<mesh_point>(cut_graph.size(), mesh_point{0, {1, 0}});
  auto dir = zero3f;
  for (auto i = 0; i < it; ++i) {
    for (auto j = 0; j < cut_graph.size(); ++j) {
      if (cut_graph[j].is_cut || !cut_graph[j].is_in_restricted)
        continue;
      auto curr_node = cut_graph[j];
      auto curr_point = mapping.at(j);
      auto relatives = find_relatives(cut_graph, j);
      if (relatives.size() == 1) {
        auto next_rel = find_relatives(cut_graph, relatives[0]);
        auto reference = dir_of_smoothing(
            mesh, topology, solver, gradients, mapping.at(next_rel[0]),
            mapping.at(relatives[0]), mapping.at(next_rel[1]), use_gradient);
        dir = dir_of_smoothing(mesh, topology, solver, dual_solver, gradients,
                               curr_point, mapping.at(relatives[0]),
                               use_gradient, reference);

      } else if (relatives.size() == 2) {
        handle_collapsing_nodes(cut_graph, relatives, mesh, mapping,
                                curr_node.vid);
        dir = dir_of_smoothing(mesh, topology, solver, gradients,
                               mapping.at(relatives[0]), curr_point,
                               mapping.at(relatives[1]), use_gradient);
      } else {

        dir = average_grad_3D(mesh, topology, dual_solver, solver, gradients,
                              relatives, mapping, curr_point);
      }

      auto magnitude = factor * length(dir) * topology.avg_edge_length;

      auto displacement = straightest_geodesic(
          solver, mesh.triangles, mesh.positions, mesh.normals,
          topology.adjacencies, topology.v2t, topology.angles,
          topology.total_angles, curr_point, normalize(dir), magnitude);

      auto new_sample = displacement.first.back();
      clean_bary(new_sample);
      subsitute[j] = new_sample;
    }
    factor *= 0.66;
    for (auto i = 0; i < subsitute.size(); ++i) {
      mapping[i] = subsitute[i];
    }
  }

  return mapping;
}

void smoothed_cut_locus(const shape_data &mesh, const shape_topology &topology,
                        const geodesic_solver &solver,
                        const dual_geodesic_solver &dual_solver,
                        vector<node> &cut_graph, const vector<vec3f> &gradients,
                        const int it, float &factor, const bool &use_gradient,
                        vector<mesh_point> &mapping) {
  auto subsitute = vector<pair<mesh_point, int>>{};
  auto dir = zero3f;
  for (auto i = 0; i < it; ++i) {
    for (auto j = 0; j < cut_graph.size(); ++j) {
      if (cut_graph[j].is_cut || !cut_graph[j].is_in_restricted)
        continue;

      auto curr_node = cut_graph[j];
      auto curr_point = mapping.at(j);
      auto relatives = find_relatives(cut_graph, j);
      if (relatives.size() == 1) {
        auto next_rel = find_relatives(cut_graph, relatives[0]);
        auto reference = dir_of_smoothing(
            mesh, topology, solver, gradients, mapping.at(next_rel[0]),
            mapping.at(relatives[0]), mapping.at(next_rel[1]), use_gradient);
        dir = dir_of_smoothing(mesh, topology, solver, dual_solver, gradients,
                               curr_point, mapping.at(relatives[0]),
                               use_gradient, reference);

      } else if (relatives.size() == 2) {
        handle_collapsing_nodes(cut_graph, relatives, mesh, mapping,
                                curr_node.vid);
        dir = dir_of_smoothing(mesh, topology, solver, gradients,
                               mapping.at(relatives[0]), curr_point,
                               mapping.at(relatives[1]), use_gradient);
      } else
        dir = average_grad_3D(mesh, topology, dual_solver, solver, gradients,
                              relatives, mapping, curr_point);

      auto magnitude = factor * length(dir) * topology.avg_edge_length;

      auto displacement = straightest_geodesic(
          solver, mesh.triangles, mesh.positions, mesh.normals,
          topology.adjacencies, topology.v2t, topology.angles,
          topology.total_angles, curr_point, normalize(dir), magnitude);
      auto new_sample = displacement.first.back();
      clean_bary(new_sample);
      subsitute.push_back({new_sample, j});
    }
    factor *= 0.66;
    for (auto i = 0; i < subsitute.size(); ++i) {
      mapping.at(subsitute[i].second) = subsitute[i].first;
    }
  }
}

vector<float> heat_geodesic(const igl::HeatGeodesicsData<double> &solver,
                            const int &source) {
  time_function();
  int V = (int)solver.Grad.cols();
  Eigen::VectorXd F(V);
  vector<float> f(V);
  Eigen::VectorXi s(1, 1);
  s << source;
  igl::heat_geodesics_solve(solver, s, F);
  for (int j = 0; j < V; ++j) {
    f[j] = F(j);
  }
  return f;
}
