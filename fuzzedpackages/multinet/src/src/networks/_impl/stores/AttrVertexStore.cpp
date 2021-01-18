/**
 * History:
 * - 2019.08.09 file created, following a restructuring of the previous library.
 */

#include "networks/_impl/stores/AttrVertexStore.hpp"

namespace uu {
namespace net {

AttrVertexStore::
AttrVertexStore()
{
    attributes_ = std::make_unique<core::AttributeStore<Vertex>>();
    attach(attributes_.get());
}


core::AttributeStore<Vertex>*
AttrVertexStore::
attr(
)
{
    return attributes_.get();
}


const core::AttributeStore<Vertex>*
AttrVertexStore::
attr(
) const
{
    return attributes_.get();
}


}
}

