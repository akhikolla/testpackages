#include "objects/Path.hpp"
#include "core/exceptions/WrongParameterException.hpp"

namespace uu {
namespace net {


Path::
Path(
    const Vertex* v0
) : super(v0)
{
}


const Vertex*
Path::
extend(
    const Edge* e
)
{
    auto new_last_vertex = super::extend(e);

    // check if the new vertex is already in the path (not the first)
    if (new_last_vertex != vertices_.front() &&
            vertex_set_.find(new_last_vertex) != vertex_set_.end())
    {
        throw core::WrongParameterException("the end-vertex of the edge is already present");
    }

    vertex_set_.insert(new_last_vertex);

    return new_last_vertex;
}

bool
Path::
is_cycle(
) const
{
    return (vertices_.front() == vertices_.back());
}
}
}

