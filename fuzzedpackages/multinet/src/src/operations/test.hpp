#ifndef UU_NET_MEASURES_TEST_H_
#define UU_NET_MEASURES_TEST_H_

#include "test.ipp"

namespace uu {
namespace net {

/**
 * Returns true is g is bipartite.
 * A bipartite graph is one for which its vertices can
 * be partitioned into two independent sets.
 * An independent set is a set of vertices V' in V where
 * for all a, b in V' a and b are not adjacent.
 * Complexity: O(n+m)
 */
template<typename G>
bool
is_bipartite(
    const G*  g
);

/**
 * Returns true is g1 and g2 are isomorphic.
 *
 * Complexity: non polynomial
 *
 * @todo NOT IMPLEMENTED
 */
template<typename G>
bool
are_isomorphic(
    const G*  g1,
    const G*  g2
);


}
}

#endif
