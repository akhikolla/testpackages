/**
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_DATASTRUCTURES_OBJECTS_NAMEDOBJECT_H_
#define UU_CORE_DATASTRUCTURES_OBJECTS_NAMEDOBJECT_H_

#include <iostream>
#include "core/datastructures/objects/Object.hpp"

namespace uu {
namespace core {

/** A generic class for objects with a name. */
class NamedObject
    : public Object
{

  public:

    NamedObject(
        const std::string& name
    );

    /** Returns a complete description of the object. */
    std::string
    to_string(
    ) const;

    /** Writes the string representation of the named object to an output stream. */
    friend std::ostream&
    operator<<(
        std::ostream&,
        const NamedObject&
    );

    /**
     * Name of the object.
     */
    const std::string name;
};

}
}

#endif
