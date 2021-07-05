#include <stdio.h>
#include <yocto/yocto_cli.h>
#include <yocto/yocto_commonio.h>
#include <yocto/yocto_geometry.h>

#include <thread>
#include <vector>
using namespace std;
#include "app.h"
#include <src/clio.h>
#include <src/cut_locus.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_opengl.h>
#include <yocto_gui/yocto_window.h>
using namespace yocto;

//
#include "editing.h"
#include "playback.h"

void set_common_uniforms(const App &app, const ogl_program *program) {
  auto &view = app.matrices.view;
  auto &projection = app.matrices.projection;
  set_uniform(program, "frame", identity4x4f);
  set_uniform(program, "view", view);
  set_uniform(program, "projection", projection);
  set_uniform(program, "eye", app.camera->frame.o);
  set_uniform(program, "envlight", (int)app.envlight);
  set_uniform(program, "gamma", app.shade_params.gamma);
  set_uniform(program, "exposure", app.shade_params.exposure);
  // set_uniform(program, "size", app.line_size);
  if (app.scene->environments.size()) {
    auto &env = app.scene->environments.front();
    if (env->envlight_diffuse)
      set_uniform(program, "envlight_irradiance", env->envlight_diffuse, 6);
    if (env->envlight_specular)
      set_uniform(program, "envlight_reflection", env->envlight_specular, 7);
    if (env->envlight_brdflut)
      set_uniform(program, "envlight_brdflut", env->envlight_brdflut, 8);
  }
}

void draw_scene(const App &app, const vec4i &viewport) {
  clear_ogl_framebuffer(vec4f{0, 0, 0, 1});

  // Draw mesh and environment.
  draw_scene(app.scene, app.camera, viewport, app.shade_params);

  if (app.show_points) {
    auto program = &app.shaders.at("points");
    bind_program(program);
    set_common_uniforms(app, program);
    set_uniform(program, "size", 3.0f * 0.0015f * app.line_size);

    set_uniform(program, "color", vec3f{0, 1, 0});

    auto draw_mesh_point = [&](const mesh_point &point, const vec3f &color) {
      if (point.face == -1)
        return;
      auto p = eval_position(app.mesh, point);
      static auto shape = ogl_shape{};
      set_points_shape(&shape, {p});
      set_uniform(program, "color", color);
      draw_shape(&shape);
    };

    for (int i = 0; i < app.eval_points.size(); i++) {
      draw_mesh_point(app.eval_points[i], {1, 1, 1});
    }
  }

  if (app.temp_levels > 0)
    draw_shape(app.temp_points[app.temp_levels]);
  auto camera_aspect = (float)viewport.z / (float)viewport.w;
  auto camera_yfov =
      camera_aspect >= 0
          ? (2 * yocto::atan(app.camera->film /
                             (camera_aspect * 2 * app.camera->lens)))
          : (2 * yocto::atan(app.camera->film / (2 * app.camera->lens)));
  auto view = frame_to_mat(inverse(app.camera->frame));
  auto projection = perspective_mat(
      camera_yfov, camera_aspect, app.shade_params.near, app.shade_params.far);

  if (app.gpu_shapes.find("edges") != app.gpu_shapes.end())
    gpu::draw_shape(app.gpu_shapes.at("edges"), app.gpu_shaders.at("points"),
                    gpu::Uniform("color", vec3f{0, 0, 0}));
  gpu::set_point_size(10);
  if (app.gpu_shapes.find("selected_points") != app.gpu_shapes.end()) {

    gpu::draw_shape(
        app.gpu_shapes.at("selected_points"), app.gpu_shaders.at("points"),
        gpu::Uniform("color", vec3f{0, 0, 1}),
        gpu::Uniform("frame", identity4x4f), gpu::Uniform("view", view),
        gpu::Uniform("projection", projection));
  }

  if (app.gpu_shapes.find("vector_field") != app.gpu_shapes.end())
    gpu::draw_shape(
        app.gpu_shapes.at("vector_field"), app.gpu_shaders.at("points"),
        gpu::Uniform("color", vec3f{0, 0, 1}),
        gpu::Uniform("frame", identity4x4f), gpu::Uniform("view", view),
        gpu::Uniform("projection", projection));

  if (app.show_edges) {
    auto program = &app.shaders.at("lines");
    bind_program(program);
    set_common_uniforms(app, program);
    set_uniform(program, "color", vec3f{0, 0, 0});
    draw_shape(&app.edges_shape);
  }

  if (app.show_branches) {
    auto program = &app.shaders.at("lines");
    bind_program(program);
    set_common_uniforms(app, program);

    set_uniform(program, "color", vec3f{0, 0, 1});
    set_uniform(program, "size", 0.001f);
    draw_shape(&app.branches_shape);
  }

  if (app.show_cobranches) {
    auto program = &app.shaders.at("lines");
    bind_program(program);
    set_common_uniforms(app, program);
    set_uniform(program, "color", vec3f{1, 0, 0});
    set_uniform(program, "size", 0.001f);
    draw_shape(&app.co_branches_shape);
  }
}

inline void sleep(int ms) {
  std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}

inline bool is_releasing(gui_button button) {
  return button.state == gui_button::state::releasing;
}
inline bool is_down(gui_button button) {
  return button.state == gui_button::state::down ||
         button.state == gui_button::state::pressing;
}

bool process_camera_move(App &app, const gui_input &input) {
  auto &camera = *app.camera;

  auto update_camera_frame = [&](frame3f &frame, float &focus, bool rotating,
                                 bool panning, bool zooming) {
    auto last_pos = input.mouse_last;
    auto mouse_pos = input.mouse_pos;
    auto mouse_left = is_down(input.mouse_left);
    auto mouse_right = is_down(input.mouse_right);
    // handle mouse and keyboard for navigation
    if (mouse_left) {
      auto dolly = 0.0f;
      auto pan = zero2f;
      auto rotate = zero2f;
      if (rotating) {
        if (mouse_left)
          rotate = (mouse_pos - last_pos) / 100.0f;
      }
      if (zooming) {
        if (mouse_right)
          dolly = (mouse_pos.y - last_pos.y) / 100.0f;
      }
      if (panning) {
        if (mouse_right)
          pan = (mouse_pos - last_pos) * focus / 200.0f;
      }
      pan.x = -pan.x;
      rotate.y = -rotate.y;
      update_turntable(frame, focus, rotate, dolly, pan);
    } else if (mouse_right) {
      auto dolly = 0.0f;
      auto pan = zero2f;
      auto rotate = zero2f;
      if (rotating) {
        if (mouse_left)
          rotate = (mouse_pos - last_pos) / 100.0f;
      }
      if (zooming) {
        if (mouse_right)
          dolly = (mouse_pos.y - last_pos.y) / 100.0f;
      }
      if (panning) {
        if (mouse_right)
          pan = (mouse_pos - last_pos) * focus / 200.0f;
      }
      pan.x = -pan.x;
      rotate.y = -rotate.y;
      update_turntable(frame, focus, rotate, dolly, pan);
    }
  };

  if (is_down(input.mouse_left)) {
    update_camera_frame(camera.frame, app.camera_focus, true, false, false);

    return true;
  }
  if (is_down(input.mouse_right)) {
    update_camera_frame(camera.frame, app.camera_focus, false, true, false);

    return true;
  }

  // Zoom-in/out by scrolling;
  float zoom = input.scroll.y * 0.1;
  if (zoom != 0) {
    update_turntable(camera.frame, app.camera_focus, zero2f, zoom, zero2f);
    return true;
  }

  return false;
}

bool process_user_input(App &app, const gui_input &input) {
  //  static bool yyy = false;

  auto mouse = input.mouse_pos;
  auto size = vec2f{(float)input.window_size.x, (float)input.window_size.y};
  mouse = vec2f{2 * (mouse.x / size.x) - 1, 1 - 2 * (mouse.y / size.y)};

  if (input.modifier_shift && is_down(input.mouse_left)) {
    // Here if pressing, but not clicked on an existing control point.
    auto point = intersect_mesh(app, mouse);
    if (point.face != -1) {
      app.control_points.push_back(point);
      if (app.added_points.size() != 0) {
        for (auto &points : app.added_points) {
          clear_shape(points->instance->shape);
        }
        app.added_points.clear();
      }
      auto pos = vector<vec3f>{eval_position(app.mesh, point)};
      update_glpoints(app, pos, "selected_points");
      return true;
    }
  }
  if (process_camera_move(app, input)) {
    update_camera_info(app, input);
    return false;
  }
  return false;
}

void update_app(App &app, const gui_input &input) {
  // process_gui_input(app, input); TODO(giacomo)

  if (is_active(app.widget))
    return;

  app.window_size = input.window_size;

  process_user_input(app, input);

  auto tasks = vector<vec2i>{};
}

void draw(const gui_input &input, void *data) {
  auto &app = *(App *)data;
  app.started = true;

  update_camera_info(app, input);

  // Do everything
  auto &t = app.playback_tick;
  if (app.playback && t < app.input_record.size()) {
    update_app(app, app.input_record[t]);
    t += 1;
  } else {
    update_app(app, input);
  }

  draw_scene(app, input.framebuffer_viewport);

  auto widget = app.widget;
  begin_widget(widget, "cutlocus");

  static vector<string> method_names = {"VTP", "Heat", "Graph"};
  draw_bullet_text(widget, "Inizialization");
  draw_combobox(widget, "Method", app.method, method_names);
  draw_intinput(widget, "Source", app.source, 0,
                (int)(app.mesh.positions.size() - 1));
  draw_checkbox(widget, "Use Harcoded Source", app.source_in_dialog);
  draw_separator(widget);
  draw_bullet_text(widget, "Cut Locus Computation");
  if (draw_button(widget, "Compute Spanning Tree")) {
    if (app.control_points.size() != 0) {
      app.source = (app.source_in_dialog)
                       ? app.source
                       : forced_vert_from_point(app.mesh.triangles,
                                                app.control_points[0]);
      app.control_points[0] = make_point_from_vert(
          app.mesh.triangles, app.topology.v2t, app.source);

    } else {
      app.source = (app.source_in_dialog) ? app.source : 0;
      app.control_points = {make_point_from_vert(app.mesh.triangles,
                                                 app.topology.v2t, app.source)};
    }
    printf("curr source is %d", app.source);

    if (app.added_points.size() != 0) {
      for (auto &point : app.added_points) {
        clear_shape(point->instance->shape);
      }
      app.added_points.clear();
    }

    add_points_shape(app, {app.mesh.positions[app.source]}, 0.003, {0, 0, 1});
    switch (app.method) {
    case VTP:
      app.field = exact_geodesic_distance(app.mesh.triangles,
                                          app.mesh.positions, app.source);
      break;
    case heat_method: {
      if (app.heat_solver.Grad.cols() > 0) {
        app.field = heat_geodesic(app.heat_solver, app.source);
      } else {
        auto [V, F] = libigl_wrapper(app.mesh.positions, app.mesh.triangles);
        igl::heat_geodesics_precompute(V, F, app.heat_solver);
        app.field = heat_geodesic(app.heat_solver, app.source);
      }
    } break;
    case graph_solver:
      app.field = compute_geodesic_distances(app.solver, {app.source});
      break;
    }

    if (app.vector_field.size() == 0) {
      app.vector_field =
          compute_grad(app.solver, app.mesh.triangles, app.mesh.positions,
                       app.mesh.normals, app.operators.Grad, app.field);
    }

    update_angle_deviation(app.mesh, app.topology, app.solver, app.vector_field,
                           app.cut_tree);
    auto min = -1;

    if (app.lap_field.size() == 0)
      std::tie(app.lap_field, min, app.lap_maximum_value) =
          compute_laplacian_with_extrema(app.operators, app.field);

    app.critical = field_maximum(app.field);
    app.lap_field[app.critical] = {flt_min, app.critical};
    spanning_tree(app.mesh, app.topology, app.solver, app.lap_field,
                  app.vector_field, app.critical, app.cut_tree);

    update_restricted(app.angle_threshold, app.lap_threshold,
                      app.use_filter_on_gradient, app.use_filter_on_laplacian,
                      app.growing, app.critical, app.cut_tree);

    if (app.thinning)
      thinning_in_restricted(app.mesh, app.solver, app.lap_field, app.cut_tree);

    if (app.clean_tree)
      cut_short_branches_in_restricted(app.cut_tree, app.depth_threshold);

    auto [branches, pos] = compute_branches(app.mesh, app.cut_tree);
    clear_shape(&app.branches_shape);
    if (pos.size() > 0) {
      set_vertex_buffer(&app.branches_shape, pos, 0);
      set_index_buffer(&app.branches_shape, branches);
    }

    app.show_branches = true;
  }
  if (draw_button(widget, "Close Loops")) {
    if (app.cut_tree.size() == 0 || app.control_points.size() == 0 ||
        app.field.size() == 0)
      return;
    app.generators.clear();

    auto genus =
        (int)(-(app.mesh.positions.size() - app.mesh.triangles.size() / 2) / 2 +
              1);
    if (genus != 0) {
      mixed_spanning_tree(app.mesh, app.topology, app.solver, app.cut_tree,
                          app.field, app.lap_field, app.lap_maximum_value,
                          app.critical, app.dist_tree, app.crossed_edges);
      spanning_cotree(app.mesh, app.topology, app.field, app.control_points[0],
                      app.cut_cotree, app.crossed_edges);

      for (auto i = 0; i < app.crossed_edges.size(); ++i) {
        for (auto j = 0; j < 3; ++j) {
          if (app.crossed_edges[i][j] == 0) {
            app.generators.push_back(vec2i{app.mesh.triangles[i][j],
                                           app.mesh.triangles[i][(j + 1) % 3]});
            auto nei = app.topology.adjacencies[i][j];
            auto k = find(app.topology.adjacencies[nei], i);
            app.crossed_edges[nei][k] = 2;
          }
        }
      }

      app.links =
          find_generators(app.mesh, app.topology, app.solver, app.dist_tree,
                          app.cut_cotree, app.generators, app.fast_loops);

      update_tree(app.links, app.lap_field, app.cut_tree);
      vector<vec3f> gen;
      for (auto points : app.generators) {
        gen.push_back(app.mesh.positions[points.x]);
        gen.push_back(app.mesh.positions[points.y]);
      }

      if (gen.size() > 0 && gen.size() < 500)
        add_points_shape(app, gen, 0.001, {0, 1, 0});
    }

    auto [branches, pos] = compute_branches(app.mesh, app.cut_tree);
    clear_shape(&app.branches_shape);
    if (pos.size() > 0) {
      set_vertex_buffer(&app.branches_shape, pos, 0);
      set_index_buffer(&app.branches_shape, branches);
    }

    app.show_branches = true;
  }
  if (draw_button(widget, "Smoothing")) {

    if (app.vector_field.size() == 0) {
      app.vector_field =
          compute_grad(app.solver, app.mesh.triangles, app.mesh.positions,
                       app.mesh.normals, app.operators.Grad, app.field);
    }

    app.mapping = smoothed_cut_locus(
        app.mesh, app.topology, app.solver, app.dual_solver, app.cut_tree,
        app.vector_field, app.it_grad_smooth, app.scaling_factor, true);
    app.scaling_factor = 1;
    smoothed_cut_locus(app.mesh, app.topology, app.solver, app.dual_solver,
                       app.cut_tree, app.vector_field, app.it_normal_smooth,
                       app.scaling_factor, false, app.mapping);

    auto [cobranches, pos] = compute_smoothed_branches(
        app.mesh, app.topology, app.cut_tree, app.mapping, true);
    if (pos.size() > 0) {
      set_vertex_buffer(&app.co_branches_shape, pos, 0);
      set_index_buffer(&app.co_branches_shape, cobranches);
    }
    app.show_cobranches = true;
    app.show_branches = false;
  }
  draw_separator(widget);
  draw_bullet_text(widget, "Filtering");
  draw_text(widget, "Policy:");
  continue_line(widget);
  if (draw_radiobutton(widget, "Growing", app.growing)) {
    app.growing = true;
    app.pruning = false;
    if (app.cut_tree.size() > 0 && app.field.size() > 0) {

      update_restricted(app.angle_threshold, app.lap_threshold,
                        app.use_filter_on_gradient, app.use_filter_on_laplacian,
                        app.growing, app.critical, app.cut_tree);
      if (app.clean_tree)
        cut_short_branches_in_restricted(app.cut_tree, app.depth_threshold);
      auto [branches, pos] = compute_branches(app.mesh, app.cut_tree);
      clear_shape(&app.branches_shape);
      if (pos.size() > 0) {
        set_vertex_buffer(&app.branches_shape, pos, 0);
        set_index_buffer(&app.branches_shape, branches);
      }
      app.show_branches = true;
    }
  }
  continue_line(widget);

  if (draw_radiobutton(widget, "Pruning", app.pruning)) {
    app.pruning = true;
    app.growing = false;
    if (app.cut_tree.size() > 0 && app.field.size() > 0) {
      update_restricted(app.angle_threshold, app.lap_threshold,
                        app.use_filter_on_gradient, app.use_filter_on_laplacian,
                        app.growing, app.critical, app.cut_tree);
      if (app.clean_tree)
        cut_short_branches_in_restricted(app.cut_tree, app.depth_threshold);

      auto [branches, pos] = compute_branches(app.mesh, app.cut_tree);
      clear_shape(&app.branches_shape);
      if (pos.size() > 0) {
        set_vertex_buffer(&app.branches_shape, pos, 0);
        set_index_buffer(&app.branches_shape, branches);
      }
      app.show_branches = true;
    }
  }
  if (draw_checkbox(widget, "Thinning", app.thinning)) {
    if (app.cut_tree.size() == 0 || app.control_points.size() == 0 ||
        app.field.size() == 0)
      return;
    if (app.thinning) {
      thinning(app.mesh, app.solver, app.lap_field, app.cut_tree);

    } else {
      for (auto i = 0; i < app.cut_tree.size(); ++i) {
        if (app.cut_tree[i].is_cut)
          app.cut_tree[i].is_cut = false;
      }
    }

    auto [branches, pos] = compute_branches(app.mesh, app.cut_tree);

    set_vertex_buffer(&app.branches_shape, pos, 0);
    set_index_buffer(&app.branches_shape, branches);
    app.show_branches = true;
  }
  if (draw_checkbox(widget, "Grad Filter", app.use_filter_on_gradient)) {
    if (app.field.size() == 0)
      return;
    if (app.cut_tree.size() == 0)
      return;
    update_restricted(app.angle_threshold, app.lap_threshold,
                      app.use_filter_on_gradient, app.use_filter_on_laplacian,
                      app.growing, app.critical, app.cut_tree);
    if (app.clean_tree)
      cut_short_branches_in_restricted(app.cut_tree, app.depth_threshold);

    auto [branches, pos] = compute_branches(app.mesh, app.cut_tree);
    clear_shape(&app.branches_shape);
    if (pos.size() > 0) {
      set_vertex_buffer(&app.branches_shape, pos, 0);
      set_index_buffer(&app.branches_shape, branches);
    }
    app.show_branches = true;
  }
  item_size(5);
  draw_floatinput(widget, "Teta_Min", app.min_angle, -pif, app.max_angle);
  continue_line(widget);
  draw_floatinput(widget, "Teta_Max", app.max_angle, app.min_angle, pif);
  item_size(15);
  if (draw_slider(widget, "Teta", app.angle_threshold, app.min_angle,
                  app.max_angle)) {
    if (app.field.size() == 0)
      return;
    if (app.cut_tree.size() == 0)
      return;
    update_restricted(app.angle_threshold, app.lap_threshold,
                      app.use_filter_on_gradient, app.use_filter_on_laplacian,
                      app.growing, app.critical, app.cut_tree);
    if (app.clean_tree)
      cut_short_branches_in_restricted(app.cut_tree, app.depth_threshold);

    auto [branches, pos] = compute_branches(app.mesh, app.cut_tree);
    clear_shape(&app.branches_shape);
    if (pos.size() > 0) {
      set_vertex_buffer(&app.branches_shape, pos, 0);
      set_index_buffer(&app.branches_shape, branches);
    }
    app.show_branches = true;
  }
  if (draw_checkbox(widget, "Lap Filter", app.use_filter_on_laplacian)) {
    if (app.field.size() == 0)
      return;
    if (app.cut_tree.size() == 0)
      return;
    update_restricted(app.angle_threshold, app.lap_threshold,
                      app.use_filter_on_gradient, app.use_filter_on_laplacian,
                      app.growing, app.critical, app.cut_tree);
    if (app.clean_tree)
      cut_short_branches_in_restricted(app.cut_tree, app.depth_threshold);

    auto [branches, pos] = compute_branches(app.mesh, app.cut_tree);
    clear_shape(&app.branches_shape);
    if (pos.size() > 0) {
      set_vertex_buffer(&app.branches_shape, pos, 0);
      set_index_buffer(&app.branches_shape, branches);
    }
    app.show_branches = true;
  }
  item_size(5);
  draw_floatinput(widget, "Lap_Min", app.min_laplacian, -1000.f,
                  app.max_laplacian);
  continue_line(widget);
  draw_floatinput(widget, "Lap_Max", app.max_laplacian, app.min_laplacian,
                  1000.f);
  item_size(15);
  if (draw_slider(widget, "Lap", app.lap_threshold, app.min_laplacian,
                  app.max_laplacian)) {
    if (app.field.size() == 0)
      return;
    if (app.cut_tree.size() == 0)
      return;
    update_restricted(app.angle_threshold, app.lap_threshold,
                      app.use_filter_on_gradient, app.use_filter_on_laplacian,
                      app.growing, app.critical, app.cut_tree);
    if (app.clean_tree)
      cut_short_branches_in_restricted(app.cut_tree, app.depth_threshold);

    auto [branches, pos] = compute_branches(app.mesh, app.cut_tree);
    clear_shape(&app.branches_shape);
    if (pos.size() > 0) {
      set_vertex_buffer(&app.branches_shape, pos, 0);
      set_index_buffer(&app.branches_shape, branches);
    }
    app.show_branches = true;
  }

  if (draw_checkbox(widget, "Cut Short Branches", app.clean_tree)) {
    if (app.clean_tree) {
      cut_short_branches_in_restricted(app.cut_tree, app.depth_threshold);
    } else {
      update_restricted(app.angle_threshold, app.lap_threshold,
                        app.use_filter_on_gradient, app.use_filter_on_laplacian,
                        app.growing, app.critical, app.cut_tree);
    }

    auto [branches, pos] = compute_branches(app.mesh, app.cut_tree);
    clear_shape(&app.branches_shape);
    set_vertex_buffer(&app.branches_shape, pos, 0);
    set_index_buffer(&app.branches_shape, branches);
    app.show_branches = true;
  }
  item_size(3);
  draw_intinput(widget, "Min", app.min_branch_length, 0, app.max_branch_length);
  continue_line(widget);
  draw_intinput(widget, "Max", app.max_branch_length, app.min_branch_length,
                100);
  item_size(15);
  if (draw_slider(widget, "Depth", app.depth_threshold, app.min_branch_length,
                  app.max_branch_length)) {

    cut_short_branches_in_restricted(app.cut_tree, app.depth_threshold);

    auto [branches, pos] = compute_branches(app.mesh, app.cut_tree);
    clear_shape(&app.branches_shape);
    if (pos.size() > 0) {
      set_vertex_buffer(&app.branches_shape, pos, 0);
      set_index_buffer(&app.branches_shape, branches);
    }
    app.show_branches = true;
  }

  draw_separator(widget);
  draw_bullet_text(widget, "Smoothing Parameters");
  draw_intinput(widget, "Gradient Iterations", app.it_grad_smooth, 0, 50);
  draw_intinput(widget, "Normal Iterations", app.it_normal_smooth, 0, 50);
  draw_separator(widget);
  draw_bullet_text(widget, "Visualization");
  item_size(15);
  if (draw_checkbox(widget, "Show Gradient", app.show_gradient)) {
    if (app.field.size() == 0)
      return;
    if (app.vector_field.size() == 0)
      app.vector_field =
          compute_grad(app.solver, app.mesh.triangles, app.mesh.positions,
                       app.mesh.normals, app.operators.Grad, app.field);
    update_glvector_field(app, app.vector_field, app.vector_size,
                          "vector_field");
  }
  if (draw_slider(widget, "vectors size", app.scale_factor, 0, 10)) {
    app.vector_size = 0.001 * app.scale_factor;
    if (app.vector_field.size() > 0)
      update_glvector_field(app, app.vector_field, app.vector_size,
                            "vector_field");
  }

  if (draw_slider(widget, "lift vector", app.lift_factor, 0, 0.05)) {
    if (app.vector_field.size() > 0 && app.lift_factor != 0)
      for (auto i = 0; i < app.vector_field.size(); ++i) {
        app.vector_field[i] += app.lift_factor * app.mesh.normals[i];
      }
    else {
      app.vector_field =
          compute_grad(app.solver, app.mesh.triangles, app.mesh.positions,
                       app.mesh.normals, app.operators.Grad, app.field);
    }
    update_glvector_field(app, app.vector_field, app.vector_size,
                          "vector_field");
  }

  draw_checkbox(widget, "show edges", app.show_edges);
  app.shade_params.faceted = app.show_edges;
  draw_separator(widget);
  if (draw_button(widget, " Reset")) {
    app.control_points = {};
    app.cut_tree = {};
    app.cut_cotree = {};
    app.dist_tree = {};
    app.generators = {};
    app.crossed_edges = {};
    app.lap_field = {};
    app.links = {};
    auto empty = vector<vec3f>{};
    app.show_branches = false;
    app.show_cobranches = false;
    clear_shape(&app.branches_shape);
    clear_shape(&app.co_branches_shape);
    app.field = {};
    app.vector_field = {};
    update_glpoints(app, empty, "selected_points2");
    update_glpoints(app, empty, "selected_points1");
    update_glpoints(app, empty, "selected_points");
    if (app.added_paths.size() != 0) {
      for (auto &path : app.added_paths) {
        clear_shape(path->instance->shape);
      }
      app.added_paths.clear();
    }
    if (app.added_points.size() != 0) {
      for (auto &points : app.added_points) {
        clear_shape(points->instance->shape);
      }

      app.added_points.clear();
    }
  }
  end_widget(widget);
}

int main(int num_args, const char *args[]) {
  auto app = App();

  string test = "";
  bool infolog = true;
  bool quiet = false;
  bool log_colors = true;
  string playback = "";
  int msaa = 1;

  auto cli = make_cli("bezier", "interactive viewer for mesh processing");
  add_argument(cli, "mesh", app.filename, "Model filenames");
  add_option(cli, "test", app.testname, "Test filename");
  // add_option(cli, "--time/--no-time", time, "Log times");
  add_option(cli, "info", infolog, "Log info");
  add_option(cli, "quiet", quiet, "Disable logs");
  add_option(cli, "colors", log_colors, "Colored logs");
  add_option(cli, "msaa", msaa, "OpenGL multisample anti-aliasing");
  add_option(cli, "playback", playback, "Playback recorded input session");
  parse_cli_and_handle_errors(cli, vector<string>{args, args + num_args});

  // Load model and init bvh for fast click-intersection.
  if (!load_mesh(app.filename, app.mesh, app.topology, app.operators,
                 app.solver, app.dual_solver, app.error))
    print_fatal(app.error);
  init_bvh(app);

  // Init window.
  auto win = new gui_window();
  win->msaa = msaa;
  init_window(win, {1080, 720}, "mesh viewer", true);
  win->user_data = &app;

  init_gpu(app, app.envlight);

  init_widget(app.widget, win);

  if (msaa > 1)
    set_ogl_msaa();

  if (playback.size()) {
    load_input_record(app.input_record, "apps/splinegui/test.recording");
    app.playback = true;
  }

  run_ui(win, draw);

  // TODO(giacomo): delete app
  clear_window(win);
}
