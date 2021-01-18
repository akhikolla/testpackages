/**
 * This file defines the attribute type to represent times.
 *
 * We have aliased `time_t` because time support might be improved in the future.
 * Currently we use the basic C++ time_t type, even if C++ support is limited.
 *
 * This file also provides basic functions to stream times, and a
 * default format to read/write/stream times.
 *
 * History:
 * - 2018.01.01 file created, adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_ATTRIBUTES_TIME_H_
#define UU_CORE_ATTRIBUTES_TIME_H_

#include <ctime>
#include <string>
#include <iostream>
#include <chrono>
#include "core/attributes/date.hpp"

namespace uu {
namespace core {

/** Time data type. */
typedef date::sys_seconds Time;


/**
 * When time values are read/written, this format is used if no other format is indicated.
 * It corresponds to YYYY-MM-DD HH:MM:SS, UTC time.
 */
const std::string kDEFAULT_TIME_FORMAT = "%Y-%m-%d %H:%M:%S %z";


/** Prints a string representation of a Time attribute value to an output stream. */
std::ostream&
operator<<(
    std::ostream& os,
    const Time& t
);


}
}

#endif
