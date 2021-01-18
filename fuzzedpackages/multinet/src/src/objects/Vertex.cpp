#include "objects/Vertex.hpp"

namespace uu {
namespace net {

Vertex::
Vertex(
    const std::string& name
) :
    NamedObject(name), key(name)
{
}

/*
std::shared_ptr<const Vertex>
Vertex::
create(
    const key_type& name
)
{
    return std::make_shared<const Vertex>(name);
}
*/

std::string
Vertex::
to_string() const
{
    return name;
}


std::ostream&
operator<<(std::ostream& os, const Vertex& v)
{
    os << v.to_string();
    return os;
}

}
}
