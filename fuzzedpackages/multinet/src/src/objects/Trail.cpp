#include "objects/Trail.hpp"
#include "core/exceptions/WrongParameterException.hpp"

namespace uu {
namespace net {


Trail::
Trail(
    const Vertex* v0
) : super(v0)
{
}


const Vertex*
Trail::
extend(
    const Edge* e
)
{
    auto new_last_vertex = super::extend(e);

    // check if the new vertex is already in the path (not the first)
    if (edge_set_.find(e) != edge_set_.end())
    {
        throw core::WrongParameterException("the edge is already present in the trail");
    }

    edge_set_.insert(e);

    return new_last_vertex;
}


}
}

