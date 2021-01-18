#include "core/attributes/conversion.hpp"
#include "core/exceptions/WrongFormatException.hpp"
#include <iomanip>
#include <iostream>
#include <vector>

namespace uu {
namespace core {

time_t
timegm(
    std::tm * timeptr
)
{
    timeptr->tm_isdst = 0;
    time_t local_time = mktime(timeptr);

    std::tm * tm_local_time_as_UTC = gmtime(&local_time);
    std::tm * tm_local_time = localtime(&local_time);

    tm_local_time->tm_isdst = 0;
    tm_local_time_as_UTC->tm_isdst = 0;
    time_t diff = mktime(tm_local_time) - mktime(tm_local_time_as_UTC);

    return local_time + diff;
}


std::string
to_string (
    const Time& v,
    const std::string& format
)
{

    time_t t = std::chrono::system_clock::to_time_t(v);
    std::tm* time = gmtime(&t);

    // std::tm time = *gmtime(&t);
    // std::ostringstream ss;
    // ss << std::put_time(&time, format.data()); ONLY AVAILABLE FROM GCC 5.0
    // return ss.str();

    char buffer [100];

    strftime (buffer, 100, format.data(), time);

    return std::string(buffer);
}


std::string
to_string (
    const Text& v
)
{
    return v.text;
}


Value<std::string>
to_string (
    const Value<Time>& v,
    const std::string& format
)
{
    std::string t = to_string(v.value, format);
    return Value<std::string>(t, v.null);
}


template <>
std::string
to_string (
    const Time& v
)
{
    return to_string(v, kDEFAULT_TIME_FORMAT);
}


template <>
Value<std::string>
to_string (
    const Value<Time>& v
)
{
    std::string t = to_string(v.value, kDEFAULT_TIME_FORMAT);
    return Value<std::string>(t, v.null);
}


double
to_double (
    const std::string& double_as_string
)
{
    std::istringstream double_val(double_as_string);
    double result;
    double_val >> result;

    if (!double_val)
    {
        throw WrongFormatException("Error converting string to double: " + double_as_string);
    }

    return result;
}


int
to_int (
    const std::string& int_as_string
)
{
    std::istringstream int_val(int_as_string);
    int result;
    int_val >> result;

    if (!int_val)
    {
        throw WrongFormatException("Error converting string to integer: " + int_as_string);
    }

    return result;
}


long
to_long (
    const std::string& long_as_string
)
{
    std::istringstream long_val(long_as_string);
    long result;
    long_val >> result;

    if (!long_val)
    {
        throw WrongFormatException("Error converting string to integer: " + long_as_string);
    }

    return result;
}


Time
epoch_to_time (
    const std::string& time_as_string
)
{
    int seconds_since_epoch = to_int(time_as_string);

    return epoch_to_time(seconds_since_epoch);
}


Time
epoch_to_time (
    int seconds_since_epoch
)
{
    Time epoch;
    std::istringstream in1{"1970-01-01 00:00:00 +0000"};
    in1 >> date::parse(kDEFAULT_TIME_FORMAT, epoch);

    std::chrono::seconds secs (seconds_since_epoch);

    return epoch + secs;

}


Time
to_time (
    const std::string& time_as_string
)
{
    return to_time(time_as_string, kDEFAULT_TIME_FORMAT);
}


Text
to_text (
    const std::string& text_as_string
)
{
    Text t;
    t.text = text_as_string;
    return t;
}


Time
to_time (
    const std::string& time_as_string,
    const std::string& format
)
{
    Time result;
    std::istringstream in{time_as_string};
    in >> date::parse(format, result);

    return result;
}


}
}
