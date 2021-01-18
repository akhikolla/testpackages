/**
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_DATASTRUCTURES_OBJECTS_OBJECT_H_
#define UU_CORE_DATASTRUCTURES_OBJECTS_OBJECT_H_

#include <iostream>

namespace uu {
namespace core {

/**
 * A generic class for objects that are part of any container defined in this
 * project, e.g., a row in a table, a vertex or an edge in a network, etc.
 */
class Object
{

  public:

    Object(
    );

    /**
     * Comparison operator: equality, based on the object reference.
     */
    bool
    operator==(
        const Object&
    ) const;

    /**
     * Comparison operator: disequality, based on the object reference.
     */
    bool
    operator!=(
        const Object&
    ) const;

    /**
     * Comparison operator: less than, bbased on the object reference.
     */
    bool
    operator<(
        const Object&
    ) const;

    /**
     * Comparison operator: higher than, based on the object reference.
     */
    bool
    operator>(
        const Object&
    ) const;

    /** Returns a string representation of the object. */
    std::string
    to_string(
    ) const;

    /** Writes the string representation of the object to an output stream. */
    friend std::ostream&
    operator<<(
        std::ostream&,
        const Object&
    );

};

}
}

#endif
