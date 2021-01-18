/**
 * This file defines the attribute types supported by the library and some basic operations on them.
 *
 * `AttributeType` lists the supported types, that are all basic C++ types. The only type
 * we have aliased is `time_t`, because time support might be improved adopting future
 * versions of C++, but also in this case we practically use the existing C(++) `time_t` type.
 *
 * This file also provides basic functions to stream attribute types.
 *
 * History:
 * - 2018.01.01 file created, adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_ATTRIBUTES_ATTRIBUTETYPE_H_
#define UU_CORE_ATTRIBUTES_ATTRIBUTETYPE_H_

#include <iostream>
#include <string>

namespace uu {
namespace core {


/** Supported attribute types. */
enum class AttributeType
{
    STRING, // std::string
    NUMERIC, // double - for back compatibility
    DOUBLE, // double
    INTEGER, // int
    TIME, // uu::core::Time
    TEXT //uu::core::Text
};


/** Returns a string representation of the input attribute type. */
std::string
to_string(
    const AttributeType& t
);


/** Prints a string representation of the attribute type to an output stream. */
std::ostream&
operator<<(
    std::ostream& os,
    const AttributeType& t
);


}
}

#endif
