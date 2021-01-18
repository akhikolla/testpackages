#include "core/attributes/Attribute.hpp"

namespace uu {
namespace core {

Attribute::
Attribute(
    const std::string& name,
    const AttributeType& type
) :
    name(name),
    type(type)
{
}

/*
std::shared_ptr<const Attribute>
Attribute::
create(
    const std::string& name,
    const AttributeType& type
)
{
    return std::make_shared<const Attribute>(name, type);
}
*/

std::unique_ptr<const Attribute>
Attribute::
create(
    const std::string& name,
    const AttributeType& type
)
{
    return std::make_unique<const Attribute>(name, type);
}

} // namespace "core"
} // namespace "uu"

