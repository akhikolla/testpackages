#include "objects/Dyad.hpp"

#include <sstream>

namespace uu {
namespace net {

Dyad::
Dyad(
    const Vertex* v1,
    const Vertex* v2
)
{
    super::insert(v1);
    super::insert(v2);
}

bool
Dyad::
operator==(
    const Dyad& comp
) const
{
    auto it1 = begin();
    auto it2 = comp.begin();

    for (size_t i = 0; i<2; i++)
    {
        if ((*it1)!=(*it2))
        {
            return false;
        }

        ++it1;
        ++it2;
    }

    return true;
}

std::set<const Vertex*>::const_iterator
Dyad::
begin() const
{
    return super::begin();
}

std::set<const Vertex*>::const_iterator
Dyad::
end() const
{
    return super::end();
}

std::set<const Vertex*>::const_iterator
Dyad::
find(const Vertex*& val) const
{
    return super::find(val);
}


std::string
Dyad::
to_string() const
{
    std::stringstream ss;

    auto it = begin();
    ss << "{" << (*it) << ",";
    it++;
    ss << (*it) << "}";
    return ss.str();
}


std::ostream&
operator<<(std::ostream& os, const Dyad& d)
{
    os << d.to_string();
    return os;
}

}
}
