#pragma once
#include <realtime/gpu.h>
#include <src/clio.h>
#include <src/cut_locus.h>
#include <unordered_set>
#include <yocto/yocto_bvh.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_shade.h>
#include <yocto_gui/yocto_window.h>

using namespace yocto;

#define PROFILE 0

inline std::unordered_set<int> make_set(size_t n) {
  auto result = std::unordered_set<int>();
  result.reserve(n);
  for (int i = 0; i < n; i++) {
    result.insert(i);
  }
  return result;
}

struct Added_Path {
  geodesic_path path;
  vec3f color;
  float radius;
  shade_instance *instance;
};

struct Added_Points {
  vector<mesh_point> points;
  vec3f color;
  float radius;
  shade_instance *instance;
};
enum struct editing_context : uint {
  is_doing_nothing,
  // is_rotating,
  is_editing_existing_curve,
  is_creating_new_curve,
  is_creating_first_curve,
};

enum prune_mode { ridge, topological, valence };

static string context_names[4] = {
    "is_doing_nothing",
    "is_editing_existing_curve",
    "is_creating_new_curve",
    "is_creating_first_curve",
};
struct Gui_Input {
  int selected_control_point = -1;
  int active_control_point = -1;
  mesh_point new_point = {};
  mesh_point first_handle = {};
  float drag_start_time = 0;
  mesh_point drag_start = {};
  mesh_point drag_current = {};
  ogl_shape *temp_points = new ogl_shape{};
  bool enforce_tangents = true;
  vec2f active_control_point_offset = {};
  int selected_spline = 0;
  float parameter_t = 0.0;
  mesh_point eval_point = {};
  mesh_point eval_point_cheap = {};
  bool translating = false;
  bool rotating = false;
  bool scaling = false;
};
struct App {
  // TEMPORARY DATA
  vector<Added_Path *> added_paths = {};
  vector<Added_Points *> added_points = {};

  bool show_branches = false;
  bool show_cobranches = false;
  bool restricted = false;
  int method = VTP;
  int initial_guess = 1;
  int prune_method = 1;
  float angle_threshold = 0.f;
  float lap_threshold = 0.f;
  float min_laplacian = -50.f;
  float max_laplacian = 50.f;
  float min_angle = 0;
  float max_angle = pif;
  int min_branch_length = 0;
  int max_branch_length = 5;
  float lap_maximum_value = flt_min;
  bool min_lap_is_inf = true;
  int depth_threshold = 3;
  int source = 0;
  bool source_in_dialog = false;
  int it_grad_smooth = 10;
  int it_normal_smooth = 10;
  bool clamped = false;
  bool use_gradient = true;
  bool clean_tree = true;
  bool thinning = true;
  bool reset_smoothing = true;
  bool use_filter_on_gradient = true;
  bool use_filter_on_laplacian = true;
  bool growing = true;
  bool pruning = false;
  bool fast_loops = true;
  bool show_gradient = false;
  std::deque<node> nodes_for_pruning = {};
  vector<vector<mesh_point>> polyline = {};
  vector<vector<int>> node_of_the_tree = {};
  vector<mesh_point> mapping = {};
  vector<float> field = {};
  float scaling_factor = 1;
  vector<float> angle_dev = {};
  vector<float> dist_from_bound = {};
  vector<vec3f> vector_field = {};
  vector<vec3f> directions = {};
  vector<pair<float, int>> lap_field = {};
  vector<vec2i> generators = {};
  vector<vector<int>> links = {};
  vector<int> cut_points = {};
  vector<vec3i> crossed_edges = {};
  vector<vec3i> dist_edge_mask = {};
  vector<int> cut_boundaries = {};
  int critical = {};
  vector<node> cut_cotree;
  vector<node> cut_tree;
  vector<node> dist_tree;
  vector<ogl_shape *> temp_points = {};
  int temp_levels = -1;
  ogl_shape *eval_point_shape = new ogl_shape{};
  Svg svg = {};
  mesh_point xxx_point = {};
  vector<mesh_point> control_points = {};
  vector<mesh_point> eval_points = {};
  vector<int> verts = {};
  int selected_vert = 0;
  igl::HeatGeodesicsData<double> heat_solver = {};
  vec2i window_size = {};
  bool started = false;
  int playback_tick = 0;
  bool playback = false;

  shape_data mesh = {};
  shape_topology topology = {};
  discr_diff_op operators = {};
  geodesic_solver solver = {};
  dual_geodesic_solver dual_solver = {};

  struct {
    mat4f view = identity4x4f;
    mat4f projection = identity4x4f;
    mat4f projection_view = identity4x4f;
  } matrices;

  bool recording_input = false;
  vector<gui_input> input_record = {};
  struct Editing_State {
    Gui_Input input = {};
  };
  editing_context context = editing_context::is_doing_nothing;
  Editing_State state = {};

  int editing_history_count = 0;
  vector<Editing_State> editing_history = {};
  int history_index = 0;

  const Gui_Input &input() const { return state.input; }
  Gui_Input &input() { return state.input; }

  void commit_state() {
    editing_history.resize(history_index + 1);
    editing_history[history_index] = state;
    history_index += 1;
    printf("%s: %d\n", "commit_state()", history_index);
  }

  string filename = "data/mesh.obj";
  string testname = "tests/test.json";
  string scene_name = "data/mesh.obj_gamma";

  float line_size = 1;
  float curve_size = 0.5;
  int scale_factor = 1;
  int type_of_scalar_field = 0;
  int type_of_strip = 0;
  float vector_size = 0.001;
  float lift_factor = 0;
  bool initial_strip = false;
  bool path_on_graph = false;
  vector<node> cut_nodes = {};
  mesh_point last_placed_point = {};

  float time_of_last_click = -1;

  float angle = 0;

  string error = "";

  bool show_edges = false;
  bool show_points = true;
  bool envlight = false;
  bool show_control_polygon = false;

  ogl_shape edges_shape = {};
  ogl_shape branches_shape = {};
  ogl_shape co_branches_shape = {};

  gui_widget *widget = new gui_widget{};
  shade_scene *scene = nullptr;
  shade_material *spline_material = nullptr;
  shade_material *mesh_material = nullptr;
  shade_shape *mesh_shape = nullptr;
  shade_camera *camera = {};
  gpu::Camera gpu_camera;
  yocto::shade_params shade_params{};
  float camera_focus;
  shape_bvh bvh = {};

  // Data stored on the gpu for rendering.
  std::unordered_map<string, ogl_shape> shapes;
  std::unordered_map<string, ogl_program> shaders;

  std::unordered_map<string, gpu::Shape> gpu_shapes;
  std::unordered_map<string, gpu::Shader> gpu_shaders;
};

#include <thread>
template <typename F> inline void parallel_for(int size, F &&f) {
  auto num_threads = min(size, 16);
  auto threads = vector<std::thread>(num_threads);
  auto batch_size = (size + num_threads - 1) / num_threads;

  auto batch = [&](int k) {
    int from = k * batch_size;
    int to = min(from + batch_size, size);
    for (int i = from; i < to; i++)
      f(i);
  };

  for (int k = 0; k < num_threads; k++) {
    threads[k] = std::thread(batch, k);
  }
  for (int k = 0; k < num_threads; k++) {
    threads[k].join();
  }
}

template <typename F> inline void serial_for(int size, F &&f) {
  for (int i = 0; i < size; i++) {
    f(i);
  }
}

void init_bvh(App &app);

shade_camera _make_framing_camera(const vector<vec3f> &positions);

void init_camera(App &app, const vec3f &from = vec3f{0, 0.5, 1.5},
                 const vec3f &to = {0, 0, 0});

vector<vec3f> make_normals(const vector<vec3i> &triangles,
                           const vector<vec3f> &positions);

void init_bvh(App &app);

ray3f camera_ray(const App &app, vec2f mouse);

vec2f screenspace_from_worldspace(App &app, const vec3f &position);

mesh_point intersect_mesh(const App &app, vec2f mouse);

void init_gpu(App &app, bool envlight);

void delete_app(App &app);

bool load_program(ogl_program *program, const string &vertex_filename,
                  const string &fragment_filename);

void set_points_shape(ogl_shape *shape, const vector<vec3f> &positions);
// void set_points_shape(ogl_shape *shape, const bezier_mesh &mesh,
//                       const vector<mesh_point> &points);
void set_mesh_shape(ogl_shape *shape, const vector<vec3i> &triangles,
                    const vector<vec3f> &positions,
                    const vector<vec3f> &normals);
void set_polyline_shape(ogl_shape *shape, const vector<vec3f> &positions);

void update_all_splines(App &app);

inline void update_camera_info(App &app, const gui_input &input) {
  auto &camera = *app.camera;
  auto viewport = input.framebuffer_viewport;
  camera.aspect = (float)viewport.z / (float)viewport.w;

  auto camera_yfov =
      (camera.aspect >= 0)
          ? (2 * yocto::atan(camera.film / (camera.aspect * 2 * camera.lens)))
          : (2 * yocto::atan(camera.film / (2 * camera.lens)));

  app.matrices.view = frame_to_mat(inverse(camera.frame));
  app.matrices.projection = perspective_mat(
      camera_yfov, camera.aspect, app.shade_params.near, app.shade_params.far);
  app.matrices.projection_view = app.matrices.projection * app.matrices.view;
}

void save_editing(const App &app, const string &filename);
bool load_editing(App &app, const string &filename);

void export_yocto_scene(const App &app, const string &name);

Added_Path *add_path_shape(App &app, const geodesic_path &path, float radius,
                           const vec3f &color);
Added_Path *add_path_shape(App &app, const vector<vec3f> &positions,
                           float radius, const vec3f &color);
Added_Points *add_points_shape(App &app, const vector<mesh_point> &points,
                               float radius, const vec3f &color);
Added_Points *add_points_shape(App &app, const vector<vec3f> &points,
                               float radius, const vec3f &color);
void update_path_shape(shade_shape *shape, const shape_data &mesh,
                       const vector<vec3i> &adjacencies,
                       const geodesic_path &path, float radius);
void update_path_shape(shade_shape *shape, const shape_data &mesh,
                       const vector<vec3f> &positions, float radius);
void update_path_shape(shade_shape *shape, const shape_data &mesh,
                       const vector<vec3f> &positions, float radius,
                       const float &treshold);
void update_points_shape(shade_shape *shape, const vector<vec3f> &positions,
                         float radius);
void update_points_shape(shade_shape *shape, const shape_data &mesh,
                         const vector<mesh_point> &points, float radius);

void update_glpoints(App &app, const vector<vec3f> &positions,
                     const string &name);

void update_glpoints(App &app, const vector<mesh_point> &points,
                     const string &name = "selected_points");

void update_glvector_field(App &app, const vector<vec3f> &vector_field,
                           float &scale, const string &name = "vector_field");

void update_glpatch(App &app, const vector<int> &faces,
                    const string &name = "patch");
