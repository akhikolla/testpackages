#include "core/datastructures/objects/NamedObject.hpp"

namespace uu {
namespace core {

NamedObject::
NamedObject(
    const std::string& name
) :
    name(name)
{

}


std::string
NamedObject::
to_string(
) const
{
    return "obj(" + name + ")";
}


std::ostream&
operator<<(
    std::ostream& os,
    const NamedObject& obj
)
{
    os << obj.to_string();
    return os;
}

}
}
