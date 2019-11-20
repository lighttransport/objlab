#pragma once

#include "params.hh"

namespace objlab {

///
/// Return true if something changed.
///
/// @param[out] texparam_changed True when texture parameter was changed.
/// @param[out] sort_pressed Alpha sorting button pressed.
///
bool alpha_window(gui_parameters *params, bool *texparam_changed, bool *sort_pressed);

bool render_window(gui_parameters *params);

}  // namespace objlab
