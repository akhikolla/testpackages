#ifndef UU_MEASURES_BETWEENNESS_H_
#define UU_MEASURES_BETWEENNESS_H_

#include <unordered_map>

namespace uu {
namespace net {

/**
 * Exact computation of betweenness centrality (Freeman, 1977) for unweighted graphs
 * using Brandes' algorithm (2001).
 */
template<typename G>
std::unordered_map<const Vertex*, double>
betweenness(
    const G* g
);

}
}

#include "betweenness.ipp" // definition

#endif
