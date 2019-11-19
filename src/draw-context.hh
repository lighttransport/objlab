#pragma once

#include "glad.h"

namespace objlab {

struct DrawObject {
  GLuint vb_id;  // vertex buffer id
  int numTriangles;
  int material_id;
};

struct DrawContext {

  std::vector<DrawObject> draw_objects;

};



} // namespace objlab

