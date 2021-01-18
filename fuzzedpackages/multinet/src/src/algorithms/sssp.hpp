#ifndef UU_NET_ALGORITHMS_SSSP_H_
#define UU_NET_ALGORITHMS_SSSP_H_

#include <unordered_set>
#include <queue>
#include "core/exceptions/assert_not_null.hpp"
#include "measures/order.hpp"
#include "objects/Vertex.hpp"
#include "objects/EdgeMode.hpp"

namespace uu {
namespace net {

/**
 * Returns the shortest path length from an input vertex to all other vertices,
 * where the length is expressed as the number of edges in the path.
 */
template<typename G>
std::vector<int>
single_source_path_length(
    const G* g,
    const Vertex* v,
    EdgeMode mode = EdgeMode::OUT
);

}
}

#import "sssp.ipp"

#endif
