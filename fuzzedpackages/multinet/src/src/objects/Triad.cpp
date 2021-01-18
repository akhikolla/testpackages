#include <sstream>
#include "objects/Triad.hpp"

namespace uu {
namespace net {

Triad::
Triad(
    const Vertex* v1,
    const Vertex* v2,
    const Vertex* v3
)
{
    super::insert(v1);
    super::insert(v2);
    super::insert(v3);
}

bool
Triad::
operator==(
    const Triad& comp
) const
{
    auto it1 = begin();
    auto it2 = comp.begin();

    for (size_t i = 0; i<3; i++)
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
Triad::
begin() const
{
    return super::begin();
}

std::set<const Vertex*>::const_iterator
Triad::
end() const
{
    return super::end();
}

std::set<const Vertex*>::const_iterator
Triad::
find(const Vertex*& val) const
{
    return super::find(val);
}


std::string
Triad::
to_string() const
{
    std::stringstream ss;

    auto it = begin();
    ss << "{" << (*it) << ",";
    it++;
    ss << (*it) << ",";
    it++;
    ss << (*it) << "}";
    return ss.str();
}


std::ostream&
operator<<(std::ostream& os, const Triad& d)
{
    os << d.to_string();
    return os;
}

}
}
