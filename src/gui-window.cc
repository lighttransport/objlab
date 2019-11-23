#include <algorithm>
#include <iostream>

#include "gui-window.hh"

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#endif

#include <glad/glad.h>

// deps/imgui
#include "imgui/imgui.h"

//#include "imgui/imgui_impl_glfw.h"
//#include "imgui/imgui_impl_opengl3.h"

//#include "ImGuizmo/ImGuizmo.h"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

namespace objlab {

bool alpha_window(gui_parameters *params, bool *texparam_changed,
                  bool *show_alpha_changed) {
  (*texparam_changed) = false;

  bool update = false;

  static const char *wrap_labels[] = {"GL_CLAMP_TO_EDGE", "GL_MIRRORED_REPEAT",
                                      "GL_REPEAT"};
  static const int wrap_value[] = {int(GL_CLAMP_TO_EDGE),
                                   int(GL_MIRRORED_REPEAT), int(GL_REPEAT)};

  const int *wrap_s_idx =
      std::find_if(wrap_value, wrap_value + 3,
                   [&params](int x) { return (x == params->texture_wrap_s); });
  const int *wrap_t_idx =
      std::find_if(wrap_value, wrap_value + 3,
                   [&params](int x) { return (x == params->texture_wrap_t); });

  int wrap_s = int(wrap_s_idx - wrap_value);
  int wrap_t = int(wrap_t_idx - wrap_value);

  if (ImGui::Begin("Alpha", &params->alpha_window_is_open)) {
    update |= ImGui::Checkbox("Enable alpha texturing",
                              &params->enable_alpha_texturing);
    ImGui::Separator();

    if (ImGui::Combo(
            "WrapS", &wrap_s, wrap_labels,
            /* num_items */ sizeof(wrap_labels) / sizeof(const char *))) {
      params->texture_wrap_s = wrap_value[wrap_s];
      (*texparam_changed) = true;
      update = true;
    }

    if (ImGui::Combo(
            "WrapT", &wrap_t, wrap_labels,
            /* num_items */ sizeof(wrap_labels) / sizeof(const char *))) {
      params->texture_wrap_t = wrap_value[wrap_t];
      (*texparam_changed) = true;
      update = true;
    }

    int radio_filter_min = (params->texture_filter_min == GL_LINEAR) ? 1 : 0;
    ImGui::Text("Filter min");
    ImGui::SameLine();
    if (ImGui::RadioButton("NEAREST##min", &radio_filter_min, 0)) {
      params->texture_filter_min = int(GL_NEAREST);
      (*texparam_changed) = true;
      update = true;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("LINEAR##min", &radio_filter_min, 1)) {
      params->texture_filter_min = int(GL_LINEAR);
      (*texparam_changed) = true;
      update = true;
    }

    int radio_filter_mag = (params->texture_filter_mag == GL_LINEAR) ? 1 : 0;
    ImGui::Text("Filter mag");
    ImGui::SameLine();
    if (ImGui::RadioButton("NEAREST##mag", &radio_filter_mag, 0)) {
      params->texture_filter_mag = int(GL_NEAREST);
      (*texparam_changed) = true;
      update = true;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("LINEAR##mag", &radio_filter_mag, 1)) {
      params->texture_filter_mag = int(GL_LINEAR);
      (*texparam_changed) = true;
      update = true;
    }

    if (ImGui::Checkbox("Show alpha", &params->texture_show_alpha)) {
      (*show_alpha_changed) = true;
      update = true;
    }
  }

  ImGui::End();

  return update;
}

bool render_window(gui_parameters *params) {
  bool update = false;

  if (ImGui::Begin("Render", &params->render_window_is_open)) {
    update |=
        ImGui::Checkbox("Enable backface cull", &params->enable_cull_face);
    update |= ImGui::Checkbox("Enable depth test", &params->enable_depth_test);
    update |= ImGui::Checkbox("Draw wireframe", &params->draw_wireframe);
    update |= ImGui::Checkbox("Show texture", &params->show_texture);
    update |= ImGui::Checkbox("Enable MSAA", &params->enable_msaa);

    ImGui::Separator();
    update |= ImGui::ColorPicker3("Background color", params->background_color);
  }

  ImGui::End();

  return update;
}

bool mesh_window(gui_parameters *params, bool *sort_pressed,
                 bool *save_pressed) {
  bool update = false;

  if (ImGui::Begin("Mesh", &params->mesh_window_is_open)) {
    ImGui::InputFloat3("view origin", params->alpha_view_origin);
    ImGui::InputFloat3("view dir", params->alpha_view_dir);
    if (ImGui::Button("Sort with above view setting.")) {
      (*sort_pressed) = true;
      update = true;
    }

    ImGui::Separator();

    if (ImGui::Button("Save sorted mesh as .obj")) {
      update = true;
      (*save_pressed) = true;
    }
  }

  ImGui::End();

  return update;
}

}  // namespace objlab
