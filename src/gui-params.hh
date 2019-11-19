#pragma once

namespace objlab {

struct gui_parameters {

  int width = 1400;
  int height = 800;

  float prevMouseX, prevMouseY;
  bool mouseLeftPressed;
  bool mouseMiddlePressed;
  bool mouseRightPressed;
  float curr_quat[4];
  float prev_quat[4];
  float eye[3], lookat[3], up[3];

};

} // namespace objlab
