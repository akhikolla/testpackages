#ifndef CLOCK_H
#define CLOCK_H

#include <type_traits>
#include <chrono>

using Clock = std::conditional<std::chrono::high_resolution_clock::is_steady,
                               std::chrono::high_resolution_clock,
                               std::chrono::steady_clock>::type;
using Seconds = std::chrono::duration<long double>;
using duration = Clock::duration;
using time_point = Clock::time_point;
using std::chrono::duration_cast;

#endif
