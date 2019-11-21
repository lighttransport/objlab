#pragma once

#include "params.hh"

namespace objlab {

///
/// @return true if something changed.
///
/// @param[out] texparam_changed True when texture parameter was changed.
/// @param[out] show_alpha_changed Show alpha was changed.
///
bool alpha_window(gui_parameters *params, bool *texparam_changed, bool *show_alpha_changed);

bool render_window(gui_parameters *params);

bool mesh_window(gui_parameters *params, bool *sort_pressed, bool *save_pressed);

}  // namespace objlab
