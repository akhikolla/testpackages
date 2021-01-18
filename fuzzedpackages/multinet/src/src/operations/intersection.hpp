#ifndef UU_OPERATIONS_INTERSECTION_H_
#define UU_OPERATIONS_INTERSECTION_H_

#include <memory>
#include <string>

namespace uu {
namespace net {

/**
 * Returns the intersection of two graphs.
 *
 * @param g1, g2 input graphs
 */
template<typename G>
std::unique_ptr<G>
graph_intersection(
    const G* g1,
    const G* g2,
    const std::string& name = ""
);


}
}

#include "intersection.ipp"

#endif
