#include "core/datastructures/objects/Object.hpp"

#include <sstream>

namespace uu {
namespace core {

Object::
Object(
)
{

}

bool
Object::
operator==(
    const Object& comp
) const
{
    return this==&comp;
}

bool
Object::
operator!=(
    const Object& comp
) const
{
    return this!=&comp;
}

bool
Object::
operator<(
    const Object& comp
) const
{
    return this<&comp;
}

bool
Object::
operator>(
    const Object& comp
) const
{
    return this>&comp;
}

std::string
Object::
to_string(
) const
{
    std::stringstream ss;
    ss << this;
    return "obj(" + ss.str() + ")";
}


std::ostream&
operator<<(
    std::ostream& os,
    const Object& obj
)
{
    os << "obj(" << (&obj) << ")";
    return os;
}

}
}
