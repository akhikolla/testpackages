#ifndef UU_OPERATIONS_SUBGRAPH_H_
#define UU_OPERATIONS_SUBGRAPH_H_

#include <memory>

namespace uu {
namespace net {

/**
 * Returns the subgraph induced by a set of vertices.
 * @param g the input graph
 * @param first, last forward iterators to the initial and final
 * positions of the sequence of const Vertex*'s.
 * The range used is [first,last).
 */
template<typename G, typename ForwardIterator>
std::unique_ptr<G>
vertex_induced_subgraph(
    const G* g,
    ForwardIterator first,
    ForwardIterator last
);


/**
 * Returns the subgraph induced by a set of edges.
 * @param g the input graph
 * @param first, last forward iterators to the initial and final
 * positions of the sequence of const Edge*'s.
 * The range used is [first,last).
 */
template<typename G, typename ForwardIterator>
std::unique_ptr<G>
edge_induced_subgraph(
    const G* g,
    ForwardIterator first,
    ForwardIterator last
);


}
}

#include "subgraph.ipp"

#endif
