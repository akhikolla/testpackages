#ifndef UU_ALGORITHMS_COMPONENTS_H_
#define UU_ALGORITHMS_COMPONENTS_H_

#include <vector>
#include <unordered_set>
#include "core/exceptions/assert_not_null.hpp"
#include "objects/Vertex.hpp"
#include "objects/EdgeMode.hpp"
#include "measures/order.hpp"
#include "algorithms/BFS.hpp"

namespace uu {
namespace net {

/**
 * Computes the components of a graph, treated as undirected.
 * @param g input graph
 * @return a vector of component ids, where the number at position i indicates the id
 * of the component of vertex g->vertices()->get_at(i)
 */
template<typename G>
std::vector<int>
components(
    const G* g
);

}
}

#import "components.ipp"

#endif
