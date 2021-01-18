#ifndef UU_MEASURES_EXPDEGREE_H_
#define UU_MEASURES_EXPDEGREE_H_

#include "objects/Vertex.hpp"
#include "objects/EdgeMode.hpp"

namespace uu {
namespace net {

/**
 * Returns the expected degree of a vertex in a probabilistic network.
 * @param g input network
 * @param v input vertex
 * @param mode to select IN, OUT, or INOUT degree
 * @return the expected degree of v in g
 */
template<typename G>
double
exp_degree(
    const G* g,
    const Vertex* v,
    const EdgeMode mode
);

}
}

#include "exp_degree.ipp"

#endif
