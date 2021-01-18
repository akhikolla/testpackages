#include "objects/Walk.hpp"

#include <sstream>
#include "objects/EdgeDir.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace net {


Walk::
Walk(
    const Vertex* v0
)
{
    core::assert_not_null(v0, "Walk", "v0");
    vertices_.push_back(v0);
}


const Vertex*
Walk::
extend(
    const Edge* e
)
{
    core::assert_not_null(e, "Walk", "e");
    auto current_last_vertex = vertices_.back();
    const Vertex* new_last_vertex;

    // check if the edge is a continuation of the walk
    if (current_last_vertex == e->v1)
    {
        new_last_vertex = e->v2;
    }

    else if (current_last_vertex == e->v2 && e->dir == EdgeDir::UNDIRECTED)
    {
        new_last_vertex = e->v1;
    }

    else
    {
        throw core::WrongParameterException("edge does not start from the last vertex in the walk");
    }

    // extend the walk
    vertices_.push_back(new_last_vertex);
    edges_.push_back(e);

    return new_last_vertex;
}


size_t
Walk::
length(
) const
{
    return edges_.size();
}



const std::list<const Vertex*>&
Walk::
vertices(
) const
{
    return vertices_;
}



const std::list<const Edge*>&
Walk::
edges(
) const
{
    return edges_;
}


std::string
Walk::
to_string() const
{
    std::stringstream ss;

    bool first = true;

    for (auto v: vertices_)
    {
        ss << (first?"":" - ") << (*v);
        first = false;
    }

    return ss.str();
}


std::ostream&
operator<<(std::ostream& os, const Walk& w)
{
    os << w.to_string();
    return os;
}

}
}

