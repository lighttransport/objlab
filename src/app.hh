#pragma once

///
/// Provides application global utility functions
///
namespace objlab {

///
/// @return Current time in seconds.
///
double GetCurrentTimeInSeconds();

///
/// Compute current frame number from elapsed time.
///
/// @param[in] start_time Start time in seconds.
/// @param[in] current_time Current time in seconds.
/// @param[in] fps Frames per second.
///
/// @return frame number
///
int ComputeCurrentFrame(double start_time, double current_time, double fps);

} // namespace objlab
