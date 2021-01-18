/**
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_ATTRIBUTES_ATTRIBUTE_H_
#define UU_CORE_ATTRIBUTES_ATTRIBUTE_H_

#include <memory>
#include <string>
#include "core/attributes/AttributeType.hpp"

namespace uu {
namespace core {

/** Meta data about an attribute. */
class Attribute :
    public std::enable_shared_from_this<Attribute>
{
  public:

    /** Name of the attribute. */
    const std::string name;

    /** Type of the attribute. */
    const AttributeType type;


    /** Constructor. */
    Attribute(
        const std::string& name,
        const AttributeType& type
    );


    /**
     * Creates a new Attribute.
     * @param name name of the attribute.
     * @param type type of the attribute.
     */
    static
    std::unique_ptr<const Attribute>
    create(
        const std::string& name,
        const AttributeType& type
    );


};

}
}

#endif
