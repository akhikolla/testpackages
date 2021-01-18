/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#include "networks/_impl/stores/EmptyEdgeStore.hpp"

namespace uu {
namespace net {


bool
EmptyEdgeStore::
is_directed(
) const
{
    return false;
}

void
EmptyEdgeStore::
erase(
    const Vertex* v
)
{
    (void)v;
}

std::string
EmptyEdgeStore::
summary(
) const
{
    return "0 edges";
}
}
}
