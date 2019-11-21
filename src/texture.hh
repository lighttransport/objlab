#pragma once

#include <vector>
#include <string>

namespace objlab {

///
/// Image contains actual image data for texturing.
///
struct Image
{
  std::string filename;
  std::vector<uint8_t> data; // Image data.
  size_t width{0};
  size_t height{0};
  int channels;
};

///
/// Texture references an array index to Image array.
/// Texture holds texture setting but does not hold actual image.
///
struct Texture
{
  uint32_t gl_tex_id{0}; // OpenGL tex id. 0 = invalid
  int image_idx{-1}; // Index to Image array.
};

} // namespace objlab
