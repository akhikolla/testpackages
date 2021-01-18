/**
 * This header defines functions to print containers and other objects.
 *
 * History:
 * - 2018.01.01 file imported from version 1.0 of the multinet library
 */

#ifndef UU_CORE_UTILS_PRETTYPRINTING_H_
#define UU_CORE_UTILS_PRETTYPRINTING_H_

#include <set>
#include <unordered_set>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

namespace uu {
namespace core {

/**
 * Converts a container to a string representation of it
 * @param first pointer to the first element in the container
 * @param last pointer to the last element in the container
 * @return a string representation of the container
 */
template <class InputIterator> std::string
set_to_string(
    InputIterator first,
    InputIterator last
)
{
    std::ostringstream ss;
    ss << "( ";

    for (; first!=last; ++first)
    {
        ss << (*first) << " ";
    }

    ss << ")";
    return ss.str();
}

/**
 * Converts a vector to a string representation of it
 * @param vec the vector
 * @return a string representation of the vector
 */
template <typename T>
std::string
to_string(
    std::vector<T> vec
)
{
    std::ostringstream ss;
    ss << "( ";

    for (T el: vec)
    {
        ss << el << " ";
    }

    ss << ")";
    return ss.str();
}

/**
 * Converts a set into a string representation of it
 * @param set the input set
 * @return a string representation of the set
 */
template <typename T>
std::string
to_string(
    std::unordered_set<T> set
)
{
    std::ostringstream ss;
    ss << "( ";

    for (T el: set)
    {
        ss << el << " ";
    }

    ss << ")";
    return ss.str();
}

/**
 * Converts a set into a string representation of it
 * @param set the input set
 * @return a string representation of the set
 */
template <typename T>
std::string
to_string(
    std::set<T> set
)
{
    std::ostringstream ss;
    ss << "( ";

    for (T el: set)
    {
        ss << el << " ";
    }

    ss << ")";
    return ss.str();
}

}
}

#endif
