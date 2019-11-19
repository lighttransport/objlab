#include "app.hh"

#if !defined(ANDROID)
// imgui must be included before glfw
#include "GLFW/glfw3.h"
#else
#include <cassert>
#endif

namespace objlab {

double GetCurrentTimeInSeconds()
{
#if defined(ANDROID)
  // TODO
  assert(0);
#else
  return glfwGetTime();
#endif
}

int ComputeCurrentFrame(double start_time, double current_time, double fps)
{
  double elapsed = current_time - start_time;

  double frame_no = elapsed * fps;

  return int(frame_no);
}

} // namespace objlab
