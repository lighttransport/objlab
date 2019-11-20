#pragma once

#include <glad/glad.h>

#include <vector>

namespace objlab {

struct DrawObject {
  GLuint vb_id;  // vertex buffer id
  int numTriangles;
  int material_id;

  std::vector<uint32_t> indices;  // TODO(LTE): index buffer
};

struct DrawContext {
  std::vector<DrawObject> draw_objects;
};

}  // namespace objlab
