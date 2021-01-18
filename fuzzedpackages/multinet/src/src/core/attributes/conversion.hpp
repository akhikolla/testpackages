/**
 * This file provides functions to convert values from/to different attribute types.
 *
 * History:
 * - 2018.01.01 file created, adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_ATTRIBUTES_CONVERSION_H_
#define UU_CORE_ATTRIBUTES_CONVERSION_H_

#include <ostream>
#include <sstream>
#include <string>
#include "core/attributes/Value.hpp"
#include "core/attributes/Text.hpp"
#include "core/attributes/Time.hpp"

namespace uu {
namespace core {

/**
 * Converts a generic type (not Time) to a string representation of it.
 * @param val value to be converted
 * @return a string representation of the input
 */
template <typename T>
std::string
to_string (
    const T& val
)
{
    std::ostringstream ss;
    ss << val;
    return ss.str();
}

/**
 * Converts a Time to a string representation of it, using the default date/time format.
 * @param v value to be converted
 * @return a string representation of the input
 */
template <>
std::string
to_string (
    const Time& v
);

/**
 * Converts a Time to a string representation of it, using a user-defined date/time format.
 * @param val value to be converted
 * @param format string format, as in the std::put_time function
 * @return a string representation of the input
 */
std::string
to_string (
    const Time& v,
    const std::string& format
);

/*
template <>
std::string
to_string (
    const Text& t
);
*/

/**
 * Converts a Value containing a generic type (not Time) to a Value containing
 * a string representation of it.
 * @param number value to be converted
 * @return a string representation of the input
 */
template <typename T>
Value<std::string>
to_string (
    const Value<T>& v
)
{
    std::ostringstream ss;
    ss << v.value;
    return Value<std::string>(ss.str(),v.null);
}


/**
 * Converts a Value containing a Time to a Value containing
 * a string representation of it.
 * @param number value to be converted
 * @return a string representation of the input
 */
template <>
Value<std::string>
to_string (
    const Value<Time>& v
);


/**
 * Converts a Value containing a Time to a Value containing
 * a string representation of it, using a user-defined date/time format.
 * @param number value to be converted
 * @param format string format, as in the std::put_time function
 * @return a string representation of the input
 */
Value<std::string>
to_string (
    const Value<Time>& v,
    const std::string& format
);


/**
 * Converts a string representation of a floating point number into its numeric value.
 * @param double_as_string a string representing a floating point number
 * @return a numerical value corresponding to the input
 * @throw WrongFormatException
 */
double
to_double(
    const std::string& double_as_string
);


/**
 * Converts a string representation of an integer into its numeric value.
 * @param int_as_string a string representing an integer
 * @return a numerical value corresponding to the input
 * @throw WrongFormatException
 */
int
to_int(
    const std::string& int_as_string
);


/**
 * Converts a string representation of a long int into its numeric value.
 * @param long_as_string a string representing an integer
 * @return a numerical value corresponding to the input
 * @throw WrongFormatException
 */
long
to_long(
    const std::string& long_as_string
);


/**
 * Converts a string representation of a time into a value of type uu::core::Time.
 * This function accepts a string representation of a long int.
 * Normally this will be interpreted as the number of seconds elapsed since 00:00
 * hours, Jan 1, 1970 UTC (i.e., a unix timestamp), but this depends on the local
 * C++ implementation. If the specific time is important and times are not just used
 * as ordinals, it is preferable to express them specifying the time format.
 *
 * @param time_as_string a string representing a time as a long int
 * @return a Time value corresponding to the input
 * @throw WrongFormatException
 */
Time
epoch_to_time(
    const std::string& seconds_since_epoch_as_string
);


/**
 * Converts a string representation of a time into a value of type uu::core::Type.
 * This function accepts a string representation of a long int.
 * Normally this will be interpreted as the number of seconds elapsed since 00:00
 * hours, Jan 1, 1970 UTC (i.e., a unix timestamp), but this depends on the local
 * C++ implementation. If the specific time is important and times are not just used
 * as ordinals, it is preferable to express them specifying the time format.
 *
 * @param time_as_string a string representing a time as a long int
 * @return a Time value corresponding to the input
 * @throw WrongFormatException
 */
Time
epoch_to_time(
    int seconds_since_epoch
);


/**
 * Converts a string representation of a time into a value of type uu::core::Type, using
 * a default format.
 *
 * @param time_as_string a string representing a time as UTC
 * @return a Time value corresponding to the input
 * @throw WrongFormatException
 */
Time
to_time(
    const std::string& time_as_string
);


/**
 * Converts a string representation of a time into a value of type uu::core::Type.
 *
 * @param time_as_string a string representing a time as UTC
 * @param the format, as described in std::get_time
 * @return a Time value corresponding to the input
 * @throw WrongFormatException
 */
Time
to_time(
    const std::string& time_as_string,
    const std::string& format
);


Text
to_text(
    const std::string& text_as_string
);


/**
 * Converts a calendar time stored in a std::tm structure
 * into a timestamp, interpreting the calendar time as UTC.
 *
 * This is (hopefully) a portable version of the inverse
 * of std::gmtime.
 */
time_t
timegm(
    std::tm * timeptr
);


}
}

#endif
