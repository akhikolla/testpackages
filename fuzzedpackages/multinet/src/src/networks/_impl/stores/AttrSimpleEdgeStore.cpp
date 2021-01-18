/**
 * History:
 * - 2019.08.09 file created, following a restructuring of the previous library.
 */

#include "networks/_impl/stores/AttrSimpleEdgeStore.hpp"

namespace uu {
namespace net {

AttrSimpleEdgeStore::
AttrSimpleEdgeStore(
    EdgeDir dir
) : SimpleEdgeStore(dir)
{
    attributes_ = std::make_unique<core::AttributeStore<Edge>>();
    attach(attributes_.get());
}


core::AttributeStore<Edge>*
AttrSimpleEdgeStore::
attr(
)
{
    return attributes_.get();
}


const core::AttributeStore<Edge>*
AttrSimpleEdgeStore::
attr(
) const
{
    return attributes_.get();
}


}
}

