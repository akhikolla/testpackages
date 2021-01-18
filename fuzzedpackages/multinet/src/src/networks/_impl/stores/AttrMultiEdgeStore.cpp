/**
 * History:
 * - 2019.08.09 file created, following a restructuring of the previous library.
 */

#include "networks/_impl/stores/AttrMultiEdgeStore.hpp"

namespace uu {
namespace net {

AttrMultiEdgeStore::
AttrMultiEdgeStore(
    EdgeDir dir
) : MultiEdgeStore(dir)
{
    attributes_ = std::make_unique<core::AttributeStore<Edge>>();
    attach(attributes_.get());
}


core::AttributeStore<Edge>*
AttrMultiEdgeStore::
attr(
)
{
    return attributes_.get();
}


const core::AttributeStore<Edge>*
AttrMultiEdgeStore::
attr(
) const
{
    return attributes_.get();
}


}
}

