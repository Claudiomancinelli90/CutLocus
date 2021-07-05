#pragma once

using namespace yocto;

inline bool is_anchor_control_point(int i) {
  if (i == -1)
    return false;
  return i % 3 == 0;
}
inline bool is_handle_control_point(int i) {
  if (i == -1)
    return false;
  return !is_anchor_control_point(i);
}

inline void set_selected_control_point(App &app, int selected,
                                       const mesh_point &pos,
                                       const vec2f &mouse) {
  auto &input = app.input();

  if (selected == -1) {
    input.selected_control_point = -1;
    input.active_control_point = -1;
    return;
  }

  auto center = screenspace_from_worldspace(app, eval_position(app.mesh, pos));

  input.active_control_point_offset = mouse - center;
  input.selected_control_point = selected;
  input.active_control_point = selected;
}

// inline mat2f parallel_transport_rotation(const bezier_mesh &mesh,
//                                          const mesh_point &start,
//                                          const mesh_point &end) {
//   auto path = compute_geodesic_path(mesh, start, end);
//   return parallel_transport_rotation(mesh.triangles, mesh.positions,
//                                      mesh.adjacencies, path);
// }
