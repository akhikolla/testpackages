/**
 * This file defines an attribute type to represent text.
 *
 * We have aliased `std::string` because text support might be improved in the future.
 *
 * This file also provides basic functions to stream text.
 *
 * History:
 * - 2018.01.03 file created
 */

#include <string>
#include <iostream>

#ifndef UU_CORE_ATTRIBUTES_TEXT_H
#define UU_CORE_ATTRIBUTES_TEXT_H

namespace uu {
namespace core {

class Text
{
  public:

    std::string text;

    /** Constructs a text object containing an empty string. */
    Text(
    ) :
        text("")
    {
    }


    /** Constructs a text object containing the input text. */
    Text(
        std::string text
    ) :
        text(text)
    {
    }
};


/** Prints a string representation of a Text attribute value to an output stream. */
std::ostream&
operator<<(
    std::ostream& os,
    const Text&
);


inline bool
operator<(const Text& lhs, const Text& rhs)
{
    return lhs.text < rhs.text;
}

inline bool
operator==(const Text& lhs, const Text& rhs)
{
    return lhs.text == rhs.text;
}

inline bool
operator!=(const Text& lhs, const Text& rhs)
{
    return lhs.text != rhs.text;
}

inline bool
operator>(const Text& lhs, const Text& rhs)
{
    return rhs < lhs;
}

inline bool
operator<=(const Text& lhs, const Text& rhs)
{
    return !(rhs < lhs);
}

inline bool
operator>=(const Text& lhs, const Text& rhs)
{
    return !(lhs < rhs);
}

}
}

#endif


