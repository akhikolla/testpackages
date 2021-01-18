/**
 * Class used to count execution times.
 *
 * History: created on 2020-04-10.
 */

#ifndef UU_CORE_UTILS_STOPWATCH_H_
#define UU_CORE_UTILS_STOPWATCH_H_

#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <random>


namespace uu {
namespace core {

class Stopwatch
{

  public:

    Stopwatch()
    {

    }

    void
    start(
    )
    {
        times.push_back(std::chrono::system_clock::now());
    }

    void
    lap(
    )
    {
        times.push_back(std::chrono::system_clock::now());
    }

    long
    millis(
        size_t lap
    )
    {
        return std::chrono::duration_cast<std::chrono::milliseconds>(times.at(lap)-times.at(lap-1)).count();
    }

    long
    sec(
        size_t lap
    )
    {
        return std::chrono::duration_cast<std::chrono::seconds>(times.at(lap)-times.at(lap-1)).count();
    }

    long
    sec_last(
    )
    {
        size_t last_lap = times.size()-1;
        return std::chrono::duration_cast<std::chrono::seconds>(times.at(last_lap)-times.at(last_lap-1)).count();
    }

  private:

    std::vector<std::chrono::time_point<std::chrono::system_clock>> times;
};

}
}

#endif
