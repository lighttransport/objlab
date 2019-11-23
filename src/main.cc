//
// Simple .obj viewer(vertex only)
//
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <chrono>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#endif

// std::filesytem in C++11
#include "ghc/filesystem.hpp"
namespace fs = ghc::filesystem;

// embeded font data for ImGui
#include "imgui/IconsIonicons.h"
#include "imgui/ionicons_embed.inc.h"
#include "imgui/roboto_mono_embed.inc.h"

// deps/imgui
#include "imgui/imgui.h"

#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include "ImGuizmo/ImGuizmo.h"

#include "glad/glad.h"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#ifdef __APPLE__
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

#include <GLFW/glfw3.h>

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#endif

#include "trackball.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#include "app.hh"
#include "draw-context.hh"
#include "face-sorter.hh"
#include "gui-window.hh"
#include "mesh.hh"
#include "params.hh"
#include "texture.hh"
#include "obj-writer.hh"

static void glfw_error_callback(int error, const char* description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

static std::string GetBaseDir(const std::string& filepath) {
  if (filepath.find_last_of("/\\") != std::string::npos)
    return filepath.substr(0, filepath.find_last_of("/\\"));
  return "";
}

static bool FileExists(const std::string& abs_filename) {
  bool ret;
  FILE* fp = fopen(abs_filename.c_str(), "rb");
  if (fp) {
    ret = true;
    fclose(fp);
  } else {
    ret = false;
  }

  return ret;
}

static void CheckGLErrors(std::string desc) {
  GLenum e = glGetError();
  if (e != GL_NO_ERROR) {
    fprintf(stderr, "OpenGL error in \"%s\": %d (%d)\n", desc.c_str(), e, e);
    exit(20);
  }
}

static void CalcNormal(float N[3], float v0[3], float v1[3], float v2[3]) {
  float v10[3];
  v10[0] = v1[0] - v0[0];
  v10[1] = v1[1] - v0[1];
  v10[2] = v1[2] - v0[2];

  float v20[3];
  v20[0] = v2[0] - v0[0];
  v20[1] = v2[1] - v0[1];
  v20[2] = v2[2] - v0[2];

  N[0] = v20[1] * v10[2] - v20[2] * v10[1];
  N[1] = v20[2] * v10[0] - v20[0] * v10[2];
  N[2] = v20[0] * v10[1] - v20[1] * v10[0];

  float len2 = N[0] * N[0] + N[1] * N[1] + N[2] * N[2];
  if (len2 > 0.0f) {
    float len = sqrtf(len2);

    N[0] /= len;
    N[1] /= len;
    N[2] /= len;
  }
}

namespace  // Local utility functions
{
struct vec3 {
  float v[3];
  vec3() {
    v[0] = 0.0f;
    v[1] = 0.0f;
    v[2] = 0.0f;
  }
};

void normalizeVector(vec3& v) {
  float len2 = v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2];
  if (len2 > 0.0f) {
    float len = sqrtf(len2);

    v.v[0] /= len;
    v.v[1] /= len;
    v.v[2] /= len;
  }
}

// Check if `mesh_t` contains smoothing group id.
bool HasSmoothingGroup(const tinyobj::shape_t& shape) {
  for (size_t i = 0; i < shape.mesh.smoothing_group_ids.size(); i++) {
    if (shape.mesh.smoothing_group_ids[i] > 0) {
      return true;
    }
  }
  return false;
}

void ComputeSmoothingNormals(const tinyobj::attrib_t& attrib,
                             const tinyobj::shape_t& shape,
                             std::map<int, vec3>& smoothVertexNormals) {
  smoothVertexNormals.clear();
  std::map<int, vec3>::iterator iter;

  for (size_t f = 0; f < shape.mesh.indices.size() / 3; f++) {
    // Get the three indexes of the face (all faces are triangular)
    tinyobj::index_t idx0 = shape.mesh.indices[3 * f + 0];
    tinyobj::index_t idx1 = shape.mesh.indices[3 * f + 1];
    tinyobj::index_t idx2 = shape.mesh.indices[3 * f + 2];

    // Get the three vertex indexes and coordinates
    int vi[3];      // indexes
    float v[3][3];  // coordinates

    for (size_t k = 0; k < 3; k++) {
      vi[0] = idx0.vertex_index;
      vi[1] = idx1.vertex_index;
      vi[2] = idx2.vertex_index;
      assert(vi[0] >= 0);
      assert(vi[1] >= 0);
      assert(vi[2] >= 0);

      v[0][k] = attrib.vertices[3 * size_t(vi[0]) + k];
      v[1][k] = attrib.vertices[3 * size_t(vi[1]) + k];
      v[2][k] = attrib.vertices[3 * size_t(vi[2]) + k];
    }

    // Compute the normal of the face
    float normal[3];
    CalcNormal(normal, v[0], v[1], v[2]);

    // Add the normal to the three vertexes
    for (size_t i = 0; i < 3; ++i) {
      iter = smoothVertexNormals.find(vi[i]);
      if (iter != smoothVertexNormals.end()) {
        // add
        iter->second.v[0] += normal[0];
        iter->second.v[1] += normal[1];
        iter->second.v[2] += normal[2];
      } else {
        smoothVertexNormals[vi[i]].v[0] = normal[0];
        smoothVertexNormals[vi[i]].v[1] = normal[1];
        smoothVertexNormals[vi[i]].v[2] = normal[2];
      }
    }

  }  // f

  // Normalize the normals, that is, make them unit vectors
  for (iter = smoothVertexNormals.begin(); iter != smoothVertexNormals.end();
       iter++) {
    normalizeVector(iter->second);
  }
}

}  // namespace

static bool LoadObjAndConvert(
    float bmin[3], float bmax[3],
    std::vector<objlab::Mesh>* meshes,                 // out
    std::vector<objlab::DrawObject>* drawObjects,      // out
    std::vector<tinyobj::material_t>& materials,       // out
    std::map<std::string, objlab::Texture>& textures,  // out
    std::vector<objlab::Image>& images,                // out
    const char* filename) {
  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t> shapes;

  meshes->clear();

  auto start_time = std::chrono::system_clock::now();

  std::string base_dir = GetBaseDir(filename);
  if (base_dir.empty()) {
    base_dir = ".";
  }
#ifdef _WIN32
  base_dir += "\\";
#else
  base_dir += "/";
#endif

  std::string warn;
  std::string err;
  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err,
                              filename, base_dir.c_str());
  if (!warn.empty()) {
    std::cout << "WARN: " << warn << std::endl;
  }
  if (!err.empty()) {
    std::cerr << err << std::endl;
  }

  auto end_time = std::chrono::system_clock::now();

  if (!ret) {
    std::cerr << "Failed to load " << filename << std::endl;
    return false;
  }

  std::chrono::duration<double, std::milli> ms = end_time - start_time;

  std::cout << "Parsing time: " << ms.count() << " [ms]\n";

  std::cout << "# of vertices  = " << (attrib.vertices.size()) / 3 << "\n";
  std::cout << "# of normals   = " << (attrib.normals.size()) / 3 << "\n";
  std::cout << "# of texcoords = " << (attrib.texcoords.size()) / 2 << "\n";
  std::cout << "# of materials = " << materials.size() << "\n";
  std::cout << "# of shapes    = " << shapes.size() << "\n";

  // Append `default` material
  materials.push_back(tinyobj::material_t());

  for (size_t i = 0; i < materials.size(); i++) {
    printf("material[%d].diffuse_texname = %s\n", int(i),
           materials[i].diffuse_texname.c_str());
    printf("material[%d].alpha_texname = %s\n", int(i),
           materials[i].alpha_texname.c_str());
  }

  // Load diffuse textures
  {
    for (size_t m = 0; m < materials.size(); m++) {
      tinyobj::material_t* mp = &materials[m];

      bool has_alpha = false;
      int alpha_w = 0, alpha_h = 0;
      std::vector<uint8_t> alpha_image;

      if (!(mp->alpha_texname.empty())) {
        int comp;

        std::string texture_filename = mp->alpha_texname;
        if (!FileExists(texture_filename)) {
          // Append base dir.
          texture_filename = base_dir + mp->alpha_texname;
          if (!FileExists(texture_filename)) {
            std::cerr << "Unable to find file: " << mp->alpha_texname
                      << std::endl;
            exit(1);
          }
        }

        unsigned char* image_data =
            stbi_load(texture_filename.c_str(), &alpha_w, &alpha_h, &comp,
                      STBI_default);
        if (!image_data) {
          std::cerr << "Unable to load texture: " << texture_filename
                    << std::endl;
          exit(1);
        }

        if (comp != 1) {
          std::cerr << "Alpha texture must be grayscale image: "
                    << texture_filename << std::endl;
          exit(1);
        }

        std::cout << "alpha_w = " << alpha_w << ", alpha_h = " << alpha_h << ", channels = " << comp << "\n";

        alpha_image.resize(size_t(alpha_w * alpha_h));
        memcpy(alpha_image.data(), image_data, size_t(alpha_w * alpha_h));

        stbi_image_free(image_data);

        has_alpha = true;
      }

      if (mp->diffuse_texname.length() > 0) {
        // Only load the texture if it is not already loaded
        if (textures.find(mp->diffuse_texname) == textures.end()) {
          GLuint texture_id;
          int w, h;
          int comp;

          std::string texture_filename = mp->diffuse_texname;
          if (!FileExists(texture_filename)) {
            // Append base dir.
            texture_filename = base_dir + mp->diffuse_texname;
            if (!FileExists(texture_filename)) {
              std::cerr << "Unable to find file: " << mp->diffuse_texname
                        << std::endl;
              exit(1);
            }
          }

          std::vector<uint8_t> pixels;
          {
            unsigned char* image_data =
                stbi_load(texture_filename.c_str(), &w, &h, &comp, STBI_default);
            if (!image_data) {
              std::cerr << "Unable to load texture: " << texture_filename
                        << std::endl;
              exit(1);
            }
            std::cout << "Loaded texture: " << texture_filename << ", w = " << w
                      << ", h = " << h << ", comp = " << comp << std::endl;

            pixels.resize(size_t(w * h * comp));
            if (comp == 4) {
              memcpy(pixels.data(), image_data, size_t(w * h * comp));

              if (has_alpha) {
                // Overwrite alpha channel with separate alpha image.
                if (alpha_w != w) {
                  std::cerr << "alpha image and color image has different image "
                               "width.\n";
                  exit(-1);
                }
                if (alpha_h != h) {
                  std::cerr << "alpha image and color image has different image "
                               "height.\n";
                  exit(-1);
                }

                for (size_t i = 0; i < size_t(w * h); i++) {
                  pixels[4 * i + 3] = alpha_image[i];
                }
              }
            } else {
              memcpy(pixels.data(), image_data, size_t(w * h * comp));
            }

            stbi_image_free(image_data);
          }

          glGenTextures(1, &texture_id);
          glBindTexture(GL_TEXTURE_2D, texture_id);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
          if (comp == 3) {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB,
                         GL_UNSIGNED_BYTE, pixels.data());
          } else if (comp == 4) {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA,
                         GL_UNSIGNED_BYTE, pixels.data());
          } else {
            assert(0);  // TODO
          }
          glBindTexture(GL_TEXTURE_2D, 0);

          objlab::Texture texture;
          texture.gl_tex_id = uint32_t(texture_id);
          texture.image_idx = int(images.size());

          // TODO(LTE): Do not create duplicated image for same filename.
          objlab::Image image;
          image.data = pixels;

          image.width = size_t(w);
          image.height = size_t(h);
          image.channels = comp;

          images.emplace_back(image);

          textures.insert(std::make_pair(mp->diffuse_texname, texture));
        }
      }
    }
  }

  bmin[0] = bmin[1] = bmin[2] = std::numeric_limits<float>::max();
  bmax[0] = bmax[1] = bmax[2] = -std::numeric_limits<float>::max();

  {
    for (size_t s = 0; s < shapes.size(); s++) {
      objlab::DrawObject o;
      std::vector<float> buffer;  // pos(3float), normal(3float), color(3float)

      objlab::Mesh mesh;

      // Check for smoothing group and compute smoothing normals
      std::map<int, vec3> smoothVertexNormals;
      if (HasSmoothingGroup(shapes[s]) > 0) {
        std::cout << "Compute smoothingNormal for shape [" << s << "]"
                  << std::endl;
        ComputeSmoothingNormals(attrib, shapes[s], smoothVertexNormals);
      }

      for (size_t f = 0; f < shapes[s].mesh.indices.size() / 3; f++) {
        tinyobj::index_t idx0 = shapes[s].mesh.indices[3 * f + 0];
        tinyobj::index_t idx1 = shapes[s].mesh.indices[3 * f + 1];
        tinyobj::index_t idx2 = shapes[s].mesh.indices[3 * f + 2];

        int current_material_id = shapes[s].mesh.material_ids[f];

        if ((current_material_id < 0) ||
            (current_material_id >= static_cast<int>(materials.size()))) {
          // Invaid material ID. Use default material.
          current_material_id =
              int(materials.size()) -
              1;  // Default material is added to the last item in `materials`.
        }
        // if (current_material_id >= materials.size()) {
        //    std::cerr << "Invalid material index: " << current_material_id <<
        //    std::endl;
        //}
        //
        float diffuse[3];
        for (size_t i = 0; i < 3; i++) {
          diffuse[i] = materials[size_t(current_material_id)].diffuse[i];
        }
        float tc[3][2];
        if (attrib.texcoords.size() > 0) {
          if ((idx0.texcoord_index < 0) || (idx1.texcoord_index < 0) ||
              (idx2.texcoord_index < 0)) {
            // face does not contain valid uv index.
            tc[0][0] = 0.0f;
            tc[0][1] = 0.0f;
            tc[1][0] = 0.0f;
            tc[1][1] = 0.0f;
            tc[2][0] = 0.0f;
            tc[2][1] = 0.0f;
          } else {
            assert(attrib.texcoords.size() >
                   size_t(2 * idx0.texcoord_index + 1));
            assert(attrib.texcoords.size() >
                   size_t(2 * idx1.texcoord_index + 1));
            assert(attrib.texcoords.size() >
                   size_t(2 * idx2.texcoord_index + 1));

            // Flip Y coord.
            tc[0][0] = attrib.texcoords[2 * size_t(idx0.texcoord_index)];
            tc[0][1] =
                1.0f - attrib.texcoords[2 * size_t(idx0.texcoord_index) + 1];
            tc[1][0] = attrib.texcoords[2 * size_t(idx1.texcoord_index)];
            tc[1][1] =
                1.0f - attrib.texcoords[2 * size_t(idx1.texcoord_index) + 1];
            tc[2][0] = attrib.texcoords[2 * size_t(idx2.texcoord_index)];
            tc[2][1] =
                1.0f - attrib.texcoords[2 * size_t(idx2.texcoord_index) + 1];
          }
        } else {
          tc[0][0] = 0.0f;
          tc[0][1] = 0.0f;
          tc[1][0] = 0.0f;
          tc[1][1] = 0.0f;
          tc[2][0] = 0.0f;
          tc[2][1] = 0.0f;
        }

        float v[3][3];
        for (size_t k = 0; k < 3; k++) {
          int f0 = idx0.vertex_index;
          int f1 = idx1.vertex_index;
          int f2 = idx2.vertex_index;
          assert(f0 >= 0);
          assert(f1 >= 0);
          assert(f2 >= 0);

          v[0][k] = attrib.vertices[3 * size_t(f0) + k];
          v[1][k] = attrib.vertices[3 * size_t(f1) + k];
          v[2][k] = attrib.vertices[3 * size_t(f2) + k];
          bmin[k] = std::min(v[0][k], bmin[k]);
          bmin[k] = std::min(v[1][k], bmin[k]);
          bmin[k] = std::min(v[2][k], bmin[k]);
          bmax[k] = std::max(v[0][k], bmax[k]);
          bmax[k] = std::max(v[1][k], bmax[k]);
          bmax[k] = std::max(v[2][k], bmax[k]);
        }

        float n[3][3];
        {
          bool invalid_normal_index = false;
          if (attrib.normals.size() > 0) {
            int nf0 = idx0.normal_index;
            int nf1 = idx1.normal_index;
            int nf2 = idx2.normal_index;

            if ((nf0 < 0) || (nf1 < 0) || (nf2 < 0)) {
              // normal index is missing from this face.
              invalid_normal_index = true;
            } else {
              for (size_t k = 0; k < 3; k++) {
                assert(size_t(3 * nf0) + k < attrib.normals.size());
                assert(size_t(3 * nf1) + k < attrib.normals.size());
                assert(size_t(3 * nf2) + k < attrib.normals.size());
                n[0][k] = attrib.normals[3 * size_t(nf0) + k];
                n[1][k] = attrib.normals[3 * size_t(nf1) + k];
                n[2][k] = attrib.normals[3 * size_t(nf2) + k];
              }
            }
          } else {
            invalid_normal_index = true;
          }

          if (invalid_normal_index && !smoothVertexNormals.empty()) {
            // Use smoothing normals
            int f0 = idx0.vertex_index;
            int f1 = idx1.vertex_index;
            int f2 = idx2.vertex_index;

            if (f0 >= 0 && f1 >= 0 && f2 >= 0) {
              n[0][0] = smoothVertexNormals[f0].v[0];
              n[0][1] = smoothVertexNormals[f0].v[1];
              n[0][2] = smoothVertexNormals[f0].v[2];

              n[1][0] = smoothVertexNormals[f1].v[0];
              n[1][1] = smoothVertexNormals[f1].v[1];
              n[1][2] = smoothVertexNormals[f1].v[2];

              n[2][0] = smoothVertexNormals[f2].v[0];
              n[2][1] = smoothVertexNormals[f2].v[1];
              n[2][2] = smoothVertexNormals[f2].v[2];

              invalid_normal_index = false;
            }
          }

          if (invalid_normal_index) {
            // compute geometric normal
            CalcNormal(n[0], v[0], v[1], v[2]);
            n[1][0] = n[0][0];
            n[1][1] = n[0][1];
            n[1][2] = n[0][2];
            n[2][0] = n[0][0];
            n[2][1] = n[0][1];
            n[2][2] = n[0][2];
          }
        }

        for (int k = 0; k < 3; k++) {
          buffer.push_back(v[k][0]);
          buffer.push_back(v[k][1]);
          buffer.push_back(v[k][2]);
          buffer.push_back(n[k][0]);
          buffer.push_back(n[k][1]);
          buffer.push_back(n[k][2]);
          // Combine normal and diffuse to get color.
          float normal_factor = 0.2f;
          float diffuse_factor = 1 - normal_factor;
          float c[3] = {n[k][0] * normal_factor + diffuse[0] * diffuse_factor,
                        n[k][1] * normal_factor + diffuse[1] * diffuse_factor,
                        n[k][2] * normal_factor + diffuse[2] * diffuse_factor};
          float len2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
          if (len2 > 0.0f) {
            float len = sqrtf(len2);

            c[0] /= len;
            c[1] /= len;
            c[2] /= len;
          }
          buffer.push_back(c[0] * 0.5f + 0.5f);
          buffer.push_back(c[1] * 0.5f + 0.5f);
          buffer.push_back(c[2] * 0.5f + 0.5f);

          buffer.push_back(tc[k][0]);
          buffer.push_back(tc[k][1]);

          mesh.vertices.push_back(v[k][0]);
          mesh.vertices.push_back(v[k][1]);
          mesh.vertices.push_back(v[k][2]);

          // facevarying normals/texcoords
          mesh.normals.push_back(n[k][0]);
          mesh.normals.push_back(n[k][1]);
          mesh.normals.push_back(n[k][2]);

          mesh.texcoords.push_back(tc[k][0]);
          mesh.texcoords.push_back(tc[k][1]);

        }

        mesh.indices.push_back(3 * uint32_t(f) + 0);
        mesh.indices.push_back(3 * uint32_t(f) + 1);
        mesh.indices.push_back(3 * uint32_t(f) + 2);

        // indices for OpenGL draw.
        o.indices.push_back(3 * uint32_t(f) + 0);
        o.indices.push_back(3 * uint32_t(f) + 1);
        o.indices.push_back(3 * uint32_t(f) + 2);

        mesh.num_verts_per_faces.push_back(3);  // triangle
      }

      o.vb_id = 0;
      o.numTriangles = 0;

      // OpenGL viewer does not support texturing with per-face material.
      if (shapes[s].mesh.material_ids.size() > 0 &&
          shapes[s].mesh.material_ids.size() > s) {
        o.material_id = shapes[s].mesh.material_ids[0];  // use the material ID
                                                         // of the first face.
      } else {
        o.material_id =
            int(materials.size()) - 1;  // = ID for default material.
      }
      printf("shape[%d] material_id %d\n", int(s), int(o.material_id));

      if (buffer.size() > 0) {
        glGenBuffers(1, &o.vb_id);
        glBindBuffer(GL_ARRAY_BUFFER, o.vb_id);
        std::cout << "vb_id = " << o.vb_id << "\n";
        glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(buffer.size() * sizeof(float)),
                     &buffer.at(0), GL_STATIC_DRAW);
        o.numTriangles = int(buffer.size() / (3 + 3 + 3 + 2) /
                             3);  // 3:vtx, 3:normal, 3:col, 2:texcoord

        printf("shape[%d] # of triangles = %d\n", static_cast<int>(s),
               o.numTriangles);
      }

      mesh.draw_object_id = int(drawObjects->size());

      drawObjects->push_back(o);
      meshes->push_back(mesh);
    }
  }

  std::cout << "bmin = " << bmin[0] << ", " << bmin[1] << ", " << bmin[2]
            << "\n";
  std::cout << "bmax = " << bmax[0] << ", " << bmax[1] << ", " << bmax[2]
            << "\n";

  return true;
}

static void ChangeTextureParameter(
    const std::map<std::string, objlab::Texture>& textures,
    const objlab::gui_parameters& params) {
  for (const auto& item : textures) {
    glBindTexture(GL_TEXTURE_2D, item.second.gl_tex_id);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, params.texture_wrap_s);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, params.texture_wrap_t);
  }

  glBindTexture(GL_TEXTURE_2D, 0);
}

static void UpdateTextures(
    const std::map<std::string, objlab::Texture>& textures,
    const std::vector<objlab::Image>& images, const bool show_alpha) {
  for (const auto& item : textures) {
    glBindTexture(GL_TEXTURE_2D, item.second.gl_tex_id);

    if (item.second.image_idx < 0) continue;  // invalid image index

    const objlab::Image& image = images[size_t(item.second.image_idx)];
    if (image.channels == 4) {
      if (show_alpha) {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, GLint(image.width),
                     GLint(image.height), 0, GL_RGBA, GL_UNSIGNED_BYTE,
                     image.data.data());
      } else {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, GLint(image.width),
                     GLint(image.height), 0, GL_RGBA, GL_UNSIGNED_BYTE,
                     image.data.data());
      }
    }
  }

  glBindTexture(GL_TEXTURE_2D, 0);
}

static void reshapeFunc(GLFWwindow* window, int w, int h) {
  std::cout << "reshape\n" << std::endl;

  int fb_w, fb_h;
  // Get actual framebuffer size.
  glfwGetFramebufferSize(window, &fb_w, &fb_h);

  glViewport(0, 0, fb_w, fb_h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0, double(w) / double(h), 0.01, 100.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  auto* param = reinterpret_cast<objlab::gui_parameters*>(
      glfwGetWindowUserPointer(window));

  param->width = w;
  param->height = h;
}

static void keyboardFunc(GLFWwindow* window, int key, int scancode, int action,
                         int mods) {
  (void)window;
  (void)scancode;
  (void)mods;

  ImGuiIO& io = ImGui::GetIO();
  if (io.WantCaptureKeyboard) {
    return;
  }

  if (action == GLFW_PRESS || action == GLFW_REPEAT) {
#if 0
    // Move camera
    float mv_x = 0, mv_y = 0, mv_z = 0;
    if (key == GLFW_KEY_K)
      mv_x += 1;
    else if (key == GLFW_KEY_J)
      mv_x += -1;
    else if (key == GLFW_KEY_L)
      mv_y += 1;
    else if (key == GLFW_KEY_H)
      mv_y += -1;
    else if (key == GLFW_KEY_P)
      mv_z += 1;
    else if (key == GLFW_KEY_N)
      mv_z += -1;
    // camera.move(mv_x * 0.05, mv_y * 0.05, mv_z * 0.05);
#endif
    // Close window
    if (key == GLFW_KEY_Q || key == GLFW_KEY_ESCAPE)
      glfwSetWindowShouldClose(window, GL_TRUE);

    // init_frame = true;
  }
}

static void clickFunc(GLFWwindow* window, int button, int action, int mods) {
  (void)mods;

  auto* params = reinterpret_cast<objlab::gui_parameters*>(
      glfwGetWindowUserPointer(window));

  if (button == GLFW_MOUSE_BUTTON_LEFT) {
    if (action == GLFW_PRESS) {
      params->mouseLeftPressed = true;
      trackball(params->prev_quat, 0.0, 0.0, 0.0, 0.0);
    } else if (action == GLFW_RELEASE) {
      params->mouseLeftPressed = false;
    }
  }
  if (button == GLFW_MOUSE_BUTTON_RIGHT) {
    if (action == GLFW_PRESS) {
      params->mouseRightPressed = true;
    } else if (action == GLFW_RELEASE) {
      params->mouseRightPressed = false;
    }
  }
  if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
    if (action == GLFW_PRESS) {
      params->mouseMiddlePressed = true;
    } else if (action == GLFW_RELEASE) {
      params->mouseMiddlePressed = false;
    }
  }
}

static void motionFunc(GLFWwindow* window, double mouse_x, double mouse_y) {
  float rotScale = 1.0f;
  float transScale = 2.0f;

  auto* params = reinterpret_cast<objlab::gui_parameters*>(
      glfwGetWindowUserPointer(window));

  if (!ImGui::GetIO().WantCaptureMouse) {
    if (params->mouseLeftPressed) {
      trackball(params->prev_quat,
                rotScale * (2.0f * params->prevMouseX - params->width) /
                    float(params->width),
                rotScale*(params->height - 2.0f * params->prevMouseY) /
                    float(params->height),
                rotScale*(2.0f * float(mouse_x) - params->width) /
                    float(params->width),
                rotScale*(params->height - 2.0f * float(mouse_y)) /
                    float(params->height));

      add_quats(params->prev_quat, params->curr_quat, params->curr_quat);
    } else if (params->mouseMiddlePressed) {
      params->eye[0] -= transScale * (float(mouse_x) - params->prevMouseX) /
                        float(params->width);
      params->lookat[0] -= transScale * (float(mouse_x) - params->prevMouseX) /
                           float(params->width);
      params->eye[1] += transScale * (float(mouse_y) - params->prevMouseY) /
                        float(params->height);
      params->lookat[1] += transScale * (float(mouse_y) - params->prevMouseY) /
                           float(params->height);
    } else if (params->mouseRightPressed) {
      params->eye[2] += transScale * (float(mouse_y) - params->prevMouseY) /
                        float(params->height);
      params->lookat[2] += transScale * (float(mouse_y) - params->prevMouseY) /
                           float(params->height);
    }
  }

  // Update mouse point
  params->prevMouseX = float(mouse_x);
  params->prevMouseY = float(mouse_y);
}

static void Draw(const std::vector<objlab::DrawObject>& drawObjects,
                 const std::vector<tinyobj::material_t>& materials,
                 const std::map<std::string, objlab::Texture>& textures,
                 bool draw_wireframe, bool show_texture) {
  glPolygonMode(GL_FRONT, GL_FILL);
  glPolygonMode(GL_BACK, GL_FILL);

  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.0, 1.0);
  GLsizei stride = (3 + 3 + 3 + 2) * sizeof(float);
  for (size_t i = 0; i < drawObjects.size(); i++) {
    objlab::DrawObject o = drawObjects[i];
    if (o.vb_id < 1) {
      continue;
    }

    glBindBuffer(GL_ARRAY_BUFFER, o.vb_id);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    if (show_texture) {
      glDisableClientState(GL_COLOR_ARRAY);
    } else {
      glEnableClientState(GL_COLOR_ARRAY);
    }
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glBindTexture(GL_TEXTURE_2D, 0);
    if (show_texture) {
      glEnable(GL_TEXTURE_2D);
      if ((o.material_id >= 0) && (o.material_id < int(materials.size()))) {
        std::string diffuse_texname =
            materials[size_t(o.material_id)].diffuse_texname;
        if (textures.find(diffuse_texname) != textures.end()) {
          glBindTexture(GL_TEXTURE_2D, textures.at(diffuse_texname).gl_tex_id);
        }
      }
    } else {
      glDisable(GL_TEXTURE_2D);
    }
    glVertexPointer(3, GL_FLOAT, stride, nullptr);
    // https://stackoverflow.com/questions/23177229/how-to-cast-int-to-const-glvoid
    glNormalPointer(
        GL_FLOAT, stride,
        reinterpret_cast<GLvoid*>(static_cast<uintptr_t>(sizeof(float) * 3)));
    glColorPointer(
        3, GL_FLOAT, stride,
        reinterpret_cast<GLvoid*>(static_cast<uintptr_t>(sizeof(float) * 6)));
    glTexCoordPointer(
        2, GL_FLOAT, stride,
        reinterpret_cast<GLvoid*>(static_cast<uintptr_t>(sizeof(float) * 9)));
    // glDrawArrays(GL_TRIANGLES, 0, 3 * o.numTriangles);
    // TODO(LTE): index bufer.
    glDrawElements(GL_TRIANGLES, GLsizei(o.indices.size()), GL_UNSIGNED_INT,
                   o.indices.data());
    CheckGLErrors("drawelements");
    glBindTexture(GL_TEXTURE_2D, 0);
  }

  if (draw_wireframe) {
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT, GL_LINE);
    glPolygonMode(GL_BACK, GL_LINE);

    glColor3f(0.0f, 0.0f, 0.4f);
    for (size_t i = 0; i < drawObjects.size(); i++) {
      objlab::DrawObject o = drawObjects[i];
      if (o.vb_id < 1) {
        continue;
      }

      glBindBuffer(GL_ARRAY_BUFFER, o.vb_id);
      glEnableClientState(GL_VERTEX_ARRAY);
      glEnableClientState(GL_NORMAL_ARRAY);
      glDisableClientState(GL_COLOR_ARRAY);
      glDisableClientState(GL_TEXTURE_COORD_ARRAY);
      glVertexPointer(3, GL_FLOAT, stride, nullptr);
      glNormalPointer(
          GL_FLOAT, stride,
          reinterpret_cast<GLvoid*>(static_cast<uintptr_t>(sizeof(float) * 3)));
      glColorPointer(
          3, GL_FLOAT, stride,
          reinterpret_cast<GLvoid*>(static_cast<uintptr_t>(sizeof(float) * 6)));
      glTexCoordPointer(
          2, GL_FLOAT, stride,
          reinterpret_cast<GLvoid*>(static_cast<uintptr_t>(sizeof(float) * 9)));

      // glDrawArrays(GL_TRIANGLES, 0, 3 * o.numTriangles);
      // TODO(LTE): index bufer.
      glDrawElements(GL_TRIANGLES, GLsizei(o.indices.size()), GL_UNSIGNED_INT,
                     o.indices.data());
      CheckGLErrors("drawarrays");
    }
  }
}

static void Init(objlab::gui_parameters* params) {
  trackball(params->curr_quat, 0, 0, 0, 0);

  params->eye[0] = 0.0f;
  params->eye[1] = 0.0f;
  params->eye[2] = 3.0f;

  params->lookat[0] = 0.0f;
  params->lookat[1] = 0.0f;
  params->lookat[2] = 0.0f;

  params->up[0] = 0.0f;
  params->up[1] = 1.0f;
  params->up[2] = 0.0f;
}

namespace {

void gui_new_frame() {
  // glfwPollEvents();
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();
}

void gl_new_frame(GLFWwindow* window, ImVec4 clear_color, bool depth_test,
                  int* display_w, int* display_h) {
  // Rendering
  glfwGetFramebufferSize(window, display_w, display_h);
  glViewport(0, 0, *display_w, *display_h);
  glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
  if (depth_test) {
    glEnable(GL_DEPTH_TEST);
  } else {
    glDisable(GL_DEPTH_TEST);
  }
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void gl_gui_end_frame(GLFWwindow* window) {
  ImGuiIO& io = ImGui::GetIO();

  // glUseProgram(0);

  ImGui::Render();

  // GLint last_program;
  // glGetIntegerv(GL_CURRENT_PROGRAM, &last_program);
  // std::cout << "last program = " << last_program << "\n";
  // glUseProgram(0);

  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

  // Update and Render additional Platform Windows
  // (Platform functions may change the current OpenGL context, so we
  // save/restore it to make it easier to paste this code elsewhere.
  //  For this specific demo app we could also call
  //  glfwMakeContextCurrent(window) directly)
  if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
    GLFWwindow* backup_current_context = glfwGetCurrentContext();
    ImGui::UpdatePlatformWindows();
    ImGui::RenderPlatformWindowsDefault();
    glfwMakeContextCurrent(backup_current_context);
  }

  // glUseProgram(last_program);

  glfwMakeContextCurrent(window);
  glfwSwapBuffers(window);
  // glFlush();

  static int fps_count = 0;
  static double currentTime = objlab::GetCurrentTimeInSeconds();
  static double previousTime = currentTime;
  static char title[256];

  fps_count++;
  currentTime = objlab::GetCurrentTimeInSeconds();
  if (currentTime - previousTime >= 1.0) {
    snprintf(title, 255, "ObjLab [%dFPS]", fps_count);
    glfwSetWindowTitle(window, title);
    fps_count = 0;
    previousTime = currentTime;
  }
}

static void initialize_imgui(GLFWwindow* window) {
  // Setup Dear ImGui context
  ImGui::CreateContext();
  auto& io = ImGui::GetIO();

  // Read .ini file from parent directory if imgui.ini does not exist in the
  // current directory.
  if (fs::exists("../imgui.ini") && !fs::exists("./imgui.ini")) {
    std::cout << "Use ../imgui.ini as Init file.\n";
    io.IniFilename = "../imgui.ini";
  }

  io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
  // io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;

  float default_font_scale = 18.0f;
  ImFontConfig roboto_config;
  strcpy(roboto_config.Name, "Roboto");
  roboto_config.SizePixels = default_font_scale;
#if 0
#if defined(__APPLE__)
  // Assume retina display. 2 is suffice
  roboto_config.OversampleH = 2;
  roboto_config.OversampleV = 2;
#else
  // 2 is a bit blurry on Windows. 4~8 gives nicer anti aliasing
  roboto_config.OversampleH = 6;
  roboto_config.OversampleV = 6;
#endif
#endif

  io.Fonts->AddFontFromMemoryCompressedTTF(roboto_mono_compressed_data,
                                           roboto_mono_compressed_size,
                                           default_font_scale, &roboto_config);

  // Load Icon fonts
  ImFontConfig ionicons_config;
  ionicons_config.MergeMode = true;
  ionicons_config.GlyphMinAdvanceX = default_font_scale;
  ionicons_config.OversampleH = 1;
  ionicons_config.OversampleV = 1;
  static const ImWchar icon_ranges[] = {ICON_MIN_II, ICON_MAX_II, 0};
  io.Fonts->AddFontFromMemoryCompressedTTF(
      ionicons_compressed_data, ionicons_compressed_size, default_font_scale,
      &ionicons_config, icon_ranges);

  ImGui::StyleColorsDark();

  // Setup Platform/Renderer bindings
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init();
}

void create_transparent_docking_area(const ImVec2 pos, const ImVec2 size,
                                     std::string name) {
  using namespace ImGui;

  const auto window_name = name + "_window";
  const auto dockspace_name = name + "_dock";

  ImGuiDockNodeFlags dockspace_flags =
      ImGuiDockNodeFlags_PassthruCentralNode |
      ImGuiDockNodeFlags_NoDockingInCentralNode;

  SetNextWindowPos(pos);
  SetNextWindowSize(size);

  ImGuiWindowFlags host_window_flags =
      ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse |
      ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
      ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoBringToFrontOnFocus |
      ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_NoBackground;

  PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
  PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
  PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));

  Begin(window_name.c_str(), nullptr, host_window_flags);
  PopStyleVar(3);  // we had 3 things added on the stack

  const ImGuiID dockspace_id = GetID(dockspace_name.c_str());
  DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), dockspace_flags);

  End();
}

void deinitialize_gui_and_window(GLFWwindow* window) {
  // Cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  glfwDestroyWindow(window);
  glfwTerminate();
}

inline void rotationY(const float angle, float* m16) {
  float c = cosf(angle);
  float s = sinf(angle);

  m16[0] = c;
  m16[1] = 0.0f;
  m16[2] = -s;
  m16[3] = 0.0f;
  m16[4] = 0.0f;
  m16[5] = 1.f;
  m16[6] = 0.0f;
  m16[7] = 0.0f;
  m16[8] = s;
  m16[9] = 0.0f;
  m16[10] = c;
  m16[11] = 0.0f;
  m16[12] = 0.f;
  m16[13] = 0.f;
  m16[14] = 0.f;
  m16[15] = 1.0f;
}

void Frustum(float left, float right, float bottom, float top, float znear,
             float zfar, float* m16) {
  float temp, temp2, temp3, temp4;
  temp = 2.0f * znear;
  temp2 = right - left;
  temp3 = top - bottom;
  temp4 = zfar - znear;
  m16[0] = temp / temp2;
  m16[1] = 0.0;
  m16[2] = 0.0;
  m16[3] = 0.0;
  m16[4] = 0.0;
  m16[5] = temp / temp3;
  m16[6] = 0.0;
  m16[7] = 0.0;
  m16[8] = (right + left) / temp2;
  m16[9] = (top + bottom) / temp3;
  m16[10] = (-zfar - znear) / temp4;
  m16[11] = -1.0f;
  m16[12] = 0.0;
  m16[13] = 0.0;
  m16[14] = (-temp * zfar) / temp4;
  m16[15] = 0.0;
}

void Perspective(float fovyInDegrees, float aspectRatio, float znear,
                 float zfar, float* m16) {
  float ymax, xmax;
  ymax = znear * std::tan(fovyInDegrees * 3.141592f / 180.0f);
  xmax = ymax * aspectRatio;
  Frustum(-xmax, xmax, -ymax, ymax, znear, zfar, m16);
}

void OrthoGraphic(const float l, float r, float b, const float t, float zn,
                  const float zf, float* m16) {
  m16[0] = 2 / (r - l);
  m16[1] = 0.0f;
  m16[2] = 0.0f;
  m16[3] = 0.0f;
  m16[4] = 0.0f;
  m16[5] = 2 / (t - b);
  m16[6] = 0.0f;
  m16[7] = 0.0f;
  m16[8] = 0.0f;
  m16[9] = 0.0f;
  m16[10] = 1.0f / (zf - zn);
  m16[11] = 0.0f;
  m16[12] = (l + r) / (l - r);
  m16[13] = (t + b) / (b - t);
  m16[14] = zn / (zn - zf);
  m16[15] = 1.0f;
}

void Cross(const float* a, const float* b, float* r) {
  r[0] = a[1] * b[2] - a[2] * b[1];
  r[1] = a[2] * b[0] - a[0] * b[2];
  r[2] = a[0] * b[1] - a[1] * b[0];
}

float Dot(const float* a, const float* b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void Normalize(const float* a, float* r) {
  float il =
      1.f / (std::sqrt(Dot(a, a)) + std::numeric_limits<float>::epsilon());
  r[0] = a[0] * il;
  r[1] = a[1] * il;
  r[2] = a[2] * il;
}

void LookAt(const float* eye, const float* at, const float* up, float* m16) {
  float X[3], Y[3], Z[3], tmp[3];

  tmp[0] = eye[0] - at[0];
  tmp[1] = eye[1] - at[1];
  tmp[2] = eye[2] - at[2];
  // Z.normalize(eye - at);
  Normalize(tmp, Z);
  Normalize(up, Y);
  // Y.normalize(up);

  Cross(Y, Z, tmp);
  // tmp.cross(Y, Z);
  Normalize(tmp, X);
  // X.normalize(tmp);

  Cross(Z, X, tmp);
  // tmp.cross(Z, X);
  Normalize(tmp, Y);
  // Y.normalize(tmp);

  m16[0] = X[0];
  m16[1] = Y[0];
  m16[2] = Z[0];
  m16[3] = 0.0f;
  m16[4] = X[1];
  m16[5] = Y[1];
  m16[6] = Z[1];
  m16[7] = 0.0f;
  m16[8] = X[2];
  m16[9] = Y[2];
  m16[10] = Z[2];
  m16[11] = 0.0f;
  m16[12] = -Dot(X, eye);
  m16[13] = -Dot(Y, eye);
  m16[14] = -Dot(Z, eye);
  m16[15] = 1.0f;
}

void EditTransform(const float* cameraView, float* cameraProjection,
                   float* matrix) {
  static ImGuizmo::OPERATION mCurrentGizmoOperation(ImGuizmo::TRANSLATE);
  static ImGuizmo::MODE mCurrentGizmoMode(ImGuizmo::LOCAL);
  static bool useSnap = false;
  static float snap[3] = {1.f, 1.f, 1.f};
  static float bounds[] = {-0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f};
  static float boundsSnap[] = {0.1f, 0.1f, 0.1f};
  static bool boundSizing = false;
  static bool boundSizingSnap = false;

  if (ImGui::IsKeyPressed(90)) mCurrentGizmoOperation = ImGuizmo::TRANSLATE;
  if (ImGui::IsKeyPressed(69)) mCurrentGizmoOperation = ImGuizmo::ROTATE;
  if (ImGui::IsKeyPressed(82))  // r Key
    mCurrentGizmoOperation = ImGuizmo::SCALE;
  if (ImGui::RadioButton("Translate",
                         mCurrentGizmoOperation == ImGuizmo::TRANSLATE))
    mCurrentGizmoOperation = ImGuizmo::TRANSLATE;
  ImGui::SameLine();
  if (ImGui::RadioButton("Rotate", mCurrentGizmoOperation == ImGuizmo::ROTATE))
    mCurrentGizmoOperation = ImGuizmo::ROTATE;
  ImGui::SameLine();
  if (ImGui::RadioButton("Scale", mCurrentGizmoOperation == ImGuizmo::SCALE))
    mCurrentGizmoOperation = ImGuizmo::SCALE;
  float matrixTranslation[3], matrixRotation[3], matrixScale[3];
  ImGuizmo::DecomposeMatrixToComponents(matrix, matrixTranslation,
                                        matrixRotation, matrixScale);
  ImGui::InputFloat3("Tr", matrixTranslation, 3);
  ImGui::InputFloat3("Rt", matrixRotation, 3);
  ImGui::InputFloat3("Sc", matrixScale, 3);
  ImGuizmo::RecomposeMatrixFromComponents(matrixTranslation, matrixRotation,
                                          matrixScale, matrix);

  if (mCurrentGizmoOperation != ImGuizmo::SCALE) {
    if (ImGui::RadioButton("Local", mCurrentGizmoMode == ImGuizmo::LOCAL))
      mCurrentGizmoMode = ImGuizmo::LOCAL;
    ImGui::SameLine();
    if (ImGui::RadioButton("World", mCurrentGizmoMode == ImGuizmo::WORLD))
      mCurrentGizmoMode = ImGuizmo::WORLD;
  }
  if (ImGui::IsKeyPressed(83)) useSnap = !useSnap;
  ImGui::Checkbox("", &useSnap);
  ImGui::SameLine();

  switch (mCurrentGizmoOperation) {
    case ImGuizmo::TRANSLATE:
      ImGui::InputFloat3("Snap", &snap[0]);
      break;
    case ImGuizmo::ROTATE:
      ImGui::InputFloat("Angle Snap", &snap[0]);
      break;
    case ImGuizmo::SCALE:
      ImGui::InputFloat("Scale Snap", &snap[0]);
      break;
    case ImGuizmo::BOUNDS:
      break;
  }
  ImGui::Checkbox("Bound Sizing", &boundSizing);
  if (boundSizing) {
    ImGui::PushID(3);
    ImGui::Checkbox("", &boundSizingSnap);
    ImGui::SameLine();
    ImGui::InputFloat3("Snap", boundsSnap);
    ImGui::PopID();
  }

  ImGuiIO& io = ImGui::GetIO();
  ImGuizmo::SetRect(0, 0, io.DisplaySize.x, io.DisplaySize.y);
  ImGuizmo::Manipulate(
      cameraView, cameraProjection, mCurrentGizmoOperation, mCurrentGizmoMode,
      matrix, nullptr, useSnap ? &snap[0] : nullptr,
      boundSizing ? bounds : nullptr, boundSizingSnap ? boundsSnap : nullptr);
}

void SortIndices(std::vector<objlab::Mesh>& meshes, // [inout]
                 std::vector<objlab::DrawObject>* draw_objects,
                 const float view_origin[3], const float view_dir[3]) {
  auto start_time = std::chrono::system_clock::now();

  for (auto& mesh : meshes) {
    std::cout << "---------------"
              << "\n";

    // Assume all triangle mesh.
    assert((mesh.indices.size() % 3) == 0);

    std::vector<uint32_t> sorted_face_indices;

    // for (size_t i = 0; i < mesh.indices.size(); i++) {
    //  std::cout << "indices = " << mesh.indices[i] << "\n";
    //}

    face_sorter::TriangleFaceCenterAccessor<float> fa(
        mesh.vertices.data(), mesh.indices.data(), mesh.indices.size() / 3);

    std::cout << "view_org = " << view_origin[0] << ", " << view_origin[1]
              << ", " << view_origin[2] << "\n";
    std::cout << "view_dir = " << view_dir[0] << ", " << view_dir[1] << ", "
              << view_dir[2] << "\n";

    face_sorter::SortByBarycentricZ<float>(mesh.indices.size() / 3, view_origin,
                                           view_dir, fa, &sorted_face_indices);

    std::cout << "sorted_face_indices = " << sorted_face_indices.size() << "\n";
    std::cout << "mesh.indices = " << mesh.indices.size() << "\n";

    int32_t draw_object_id = mesh.draw_object_id;
    assert(draw_object_id >= 0);

    // reorder vetex indices.
    std::vector<uint32_t> sorted_indices;
    for (size_t i = 0; i < sorted_face_indices.size(); i++) {
      size_t idx = sorted_face_indices[i];

      uint32_t f0 = mesh.indices[3 * idx + 0];
      uint32_t f1 = mesh.indices[3 * idx + 1];
      uint32_t f2 = mesh.indices[3 * idx + 2];

      sorted_indices.push_back(f0);
      sorted_indices.push_back(f1);
      sorted_indices.push_back(f2);
      // std::cout << "sort = " << f0 << ", " << f1 << ", " << f2 << "\n";
    }

    objlab::DrawObject& o = (*draw_objects)[size_t(draw_object_id)];

    assert(sorted_indices.size() == o.indices.size());

    // update indices
    o.indices = sorted_indices;

    // also store sorted indices to Mesh
    mesh.sorted_indices = sorted_indices;

  }

  auto end_time = std::chrono::system_clock::now();

  std::chrono::duration<double, std::milli> ms = end_time - start_time;

  std::cout << "Sort time: " << ms.count() << " [ms]\n";
}

}  // namespace

int main(int argc, char** argv) {
#if defined(_MSC_VER)
	std::string obj_filename = "../../../models/cornell_box.obj";
#else
  std::string obj_filename = "../models/cornell_box.obj";
#endif

  if (argc < 2) {
    std::cout << "Needs input.obj\n" << std::endl;
    std::cout << "  Use default value: " << obj_filename << "\n";
  }

  if (argc > 1) {
    obj_filename = argv[1];
  }

  objlab::gui_parameters gui_params;

  Init(&gui_params);

  glfwSetErrorCallback(glfw_error_callback);

  if (!glfwInit()) {
    std::cerr << "Failed to initialize GLFW." << std::endl;
    return -1;
  }

  // MSAA
  glfwWindowHint(GLFW_SAMPLES, 16);

#if 0
    // Decide GL+GLSL versions
#if defined(__APPLE__)
    // GL 3.2 + GLSL 150
    const char* glsl_version = "#version 150";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac
#else
    // GL 3.0 + GLSL 130
    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
#endif
#endif

  GLFWwindow* window = glfwCreateWindow(gui_params.width, gui_params.height,
                                        "Obj viewer", nullptr, nullptr);
  if (window == nullptr) {
    std::cerr << "Failed to open GLFW window. " << std::endl;
    glfwTerminate();
    return 1;
  }

  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);

  glfwSetWindowUserPointer(window, &gui_params);

  // Callback
  glfwSetWindowSizeCallback(window, reshapeFunc);
  glfwSetKeyCallback(window, keyboardFunc);
  glfwSetMouseButtonCallback(window, clickFunc);
  glfwSetCursorPosCallback(window, motionFunc);

#if 1
  if (gladLoadGL() == 0) {
    std::cerr << "Failed to initialize GLAD.\n";
    return -1;
  }
#else
  if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress))) {
    std::cerr << "Failed to initialize GLAD.\n";
    return -1;
  }
#endif

  if (!GLAD_GL_VERSION_3_0) {
    std::cerr << "OpenGL 3.0 context not available.\n";
    return -1;
  }

  initialize_imgui(window);

  reshapeFunc(window, gui_params.width, gui_params.height);

  objlab::DrawContext draw_ctx;

  float bmin[3], bmax[3];
  std::vector<tinyobj::material_t> materials;
  std::map<std::string, objlab::Texture> textures;
  std::vector<objlab::Image> images;
  std::vector<objlab::Mesh> meshes;

  if (false == LoadObjAndConvert(bmin, bmax, &meshes, &draw_ctx.draw_objects,
                                 materials, textures, images,
                                 obj_filename.c_str())) {
    return -1;
  }

  float maxExtent = 0.5f * (bmax[0] - bmin[0]);
  if (maxExtent < 0.5f * (bmax[1] - bmin[1])) {
    maxExtent = 0.5f * (bmax[1] - bmin[1]);
  }
  if (maxExtent < 0.5f * (bmax[2] - bmin[2])) {
    maxExtent = 0.5f * (bmax[2] - bmin[2]);
  }

  float objectMatrix[16] = {1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f,
                            0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f};

  static const float identityMatrix[16] = {1.f, 0.f, 0.f, 0.f, 0.f, 1.f,
                                           0.f, 0.f, 0.f, 0.f, 1.f, 0.f,
                                           0.f, 0.f, 0.f, 1.f};

  float cameraView[16] = {1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f,
                          0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f};

  float cameraProjection[16];

  // Camera projection
  bool isPerspective = false;
  float fov = 27.f;
  float viewWidth = 10.f;  // for orthographic
  float camYAngle = 165.f / 180.f * 3.14159f;
  float camXAngle = 0.f / 180.f * 3.14159f;
  float camDistance = 8.f;
  rotationY(0.f, objectMatrix);

  bool firstFrame = true;

  while (glfwWindowShouldClose(window) == GL_FALSE) {
    int display_w, display_h;

    glfwPollEvents();

    gui_new_frame();

#if 0
    ImGuiIO& io = ImGui::GetIO();
    if (isPerspective) {
      Perspective(fov, io.DisplaySize.x / io.DisplaySize.y, 0.1f, 100.f,
                  cameraProjection);
    } else {
      float viewHeight = viewWidth * io.DisplaySize.y / io.DisplaySize.x;
      OrthoGraphic(-viewWidth, viewWidth, -viewHeight, viewHeight, -viewWidth,
                   viewWidth, cameraProjection);
    }
    ImGuizmo::SetOrthographic(!isPerspective);

    ImGuizmo::BeginFrame();

    ImGui::SetNextWindowPos(ImVec2(1024, 100));
    ImGui::SetNextWindowSize(ImVec2(256, 256));

    // create a window and insert the inspector
    ImGui::SetNextWindowPos(ImVec2(10, 10));
    ImGui::SetNextWindowSize(ImVec2(320, 340));
    ImGui::Begin("Editor");
    ImGui::Text("Camera");
    bool viewDirty = false;
    if (ImGui::RadioButton("Perspective", isPerspective)) isPerspective = true;
    ImGui::SameLine();
    if (ImGui::RadioButton("Orthographic", !isPerspective))
      isPerspective = false;
    if (isPerspective) {
      ImGui::SliderFloat("Fov", &fov, 20.f, 110.f);
    } else {
      ImGui::SliderFloat("Ortho width", &viewWidth, 1, 20);
    }
    viewDirty |= ImGui::SliderAngle("Camera X", &camXAngle, 0.f, 179.f);
    viewDirty |= ImGui::SliderAngle("Camera Y", &camYAngle);
    viewDirty |= ImGui::SliderFloat("Distance", &camDistance, 1.f, 10.f);

    if (viewDirty || firstFrame) {
      float eye[] = {cosf(camYAngle) * cosf(camXAngle) * camDistance + 2.f,
                     sinf(camXAngle) * camDistance,
                     sinf(camYAngle) * cosf(camXAngle) * camDistance};
      float at[] = {2.f, 0.f, 0.f};
      float up[] = {0.f, 1.f, 0.f};
      LookAt(eye, at, up, cameraView);
      firstFrame = false;
    }
    //ImGuizmo::DrawCube(cameraView, cameraProjection, objectMatrix);
    //ImGuizmo::DrawGrid(cameraView, cameraProjection, identityMatrix, 10.f);

    ImGui::Text("X: %f Y: %f", double(io.MousePos.x), double(io.MousePos.y));
    ImGui::Separator();
    EditTransform(cameraView, cameraProjection, objectMatrix);
    ImGui::End(); // 'Editor'

    {
      // Create docking area
      const auto display_size = ImGui::GetIO().DisplaySize;
      // Use (0, 20) if you need menu bar
      create_transparent_docking_area(ImVec2(0, 0),
                                      ImVec2(display_size.x, display_size.y),
                                      "main_dockspace");
    }

    ImGuizmo::ViewManipulate(cameraView, camDistance,
                             ImVec2(io.DisplaySize.x - 128, 0),
                             ImVec2(128, 128), 0x10101010);
#endif

    {
      bool sort_pressed = false;
      bool save_pressed = false;
      if (mesh_window(&gui_params, &sort_pressed, &save_pressed)) {

        if (sort_pressed) {
          SortIndices(meshes, &draw_ctx.draw_objects,
                      gui_params.alpha_view_origin, gui_params.alpha_view_dir);
        }


        if (save_pressed) {
          if (!SaveMeshAsObj(meshes, gui_params.output_obj_basename)) {
            std::cerr << "Failed to save .obj.\n";
          } else {
            std::cout << "Saved .obj : " << gui_params.output_obj_basename << ".obj\n";
          }
        }

      }
    }

    {
      bool texparam_changed = false;
      bool show_alpha_changed = false;
      if (alpha_window(&gui_params, &texparam_changed,
                       &show_alpha_changed)) {
        if (texparam_changed) {
          ChangeTextureParameter(textures, gui_params);
        }

        if (show_alpha_changed) {
          UpdateTextures(textures, images, gui_params.texture_show_alpha);
        }
      }
    }

    if (render_window(&gui_params)) {
    }

    ImVec4 background_color = {gui_params.background_color[0],
                               gui_params.background_color[1],
                               gui_params.background_color[2], 1.0f};
    gl_new_frame(window, background_color, gui_params.enable_depth_test,
                 &display_w, &display_h);

    {
      if (gui_params.enable_msaa) {
        glEnable(GL_MULTISAMPLE);
      } else {
        glDisable(GL_MULTISAMPLE);
      }

      glEnable(GL_TEXTURE_2D);

      if (gui_params.enable_alpha_texturing) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
      } else {
        glDisable(GL_BLEND);

        if (gui_params.enable_depth_test) {
          glEnable(GL_DEPTH_TEST);
        } else {
          glDisable(GL_DEPTH_TEST);
        }

        if (gui_params.enable_cull_face) {
          glEnable(GL_CULL_FACE);
        } else {
          glDisable(GL_CULL_FACE);
        }
      }

#if 1
      // camera & rotate
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      GLfloat mat[4][4];
      gluLookAt(double(gui_params.eye[0]), double(gui_params.eye[1]),
                double(gui_params.eye[2]), double(gui_params.lookat[0]),
                double(gui_params.lookat[1]), double(gui_params.lookat[2]),
                double(gui_params.up[0]), double(gui_params.up[1]),
                double(gui_params.up[2]));
      build_rotmatrix(mat, gui_params.curr_quat);
      glMultMatrixf(&mat[0][0]);
#endif

      // Fit to -1, 1
      glScalef(1.0f / maxExtent, 1.0f / maxExtent, 1.0f / maxExtent);

      // Centerize object.
      glTranslatef(-0.5f * (bmax[0] + bmin[0]), -0.5f * (bmax[1] + bmin[1]),
                   -0.5f * (bmax[2] + bmin[2]));

      Draw(draw_ctx.draw_objects, materials, textures,
           gui_params.draw_wireframe, gui_params.show_texture);

      glEnable(GL_DEPTH_TEST);
      glEnable(GL_CULL_FACE);
      glDisable(GL_TEXTURE_2D);
    }

    gl_gui_end_frame(window);
  }

  deinitialize_gui_and_window(window);

  return EXIT_SUCCESS;
}
