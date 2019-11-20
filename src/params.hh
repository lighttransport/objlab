#pragma once

#include <glad/glad.h>

namespace objlab {

struct gui_parameters {

  int width = 1400;
  int height = 800;

  float prevMouseX = 0.0f, prevMouseY = 0.0f;
  bool mouseLeftPressed = false;
  bool mouseMiddlePressed = false;
  bool mouseRightPressed = false;
  float curr_quat[4] = {0.0f, 0.0f, 0.0f, 1.0f};
  float prev_quat[4] = {0.0f, 0.0f, 0.0f, 1.0f};
  float eye[3], lookat[3], up[3];

  // RGB
  float background_color[3] = {0.1f, 0.15f, 0.2f};

  bool alpha_window_is_open = true;
  bool enable_alpha_texturing = true;

  bool render_window_is_open = true;
  bool enable_cull_face = true;
  bool enable_depth_test = true;
  bool draw_wireframe = true;

  // polygon sorting for alpha
  // (coordinate is defined in right-handed)
  float alpha_view_origin[3] = {0.0f, 0.0f, 1.0f};
  float alpha_view_dir[3] = {0.0f, 0.0f, -1.0f};

  // texture parameters
  int texture_wrap_s = GL_CLAMP_TO_EDGE;
  int texture_wrap_t = GL_CLAMP_TO_EDGE;
};

} // namespace objlab
