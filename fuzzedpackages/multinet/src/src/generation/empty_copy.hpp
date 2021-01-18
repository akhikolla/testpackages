#ifndef UU_OPERATIONS_EMPTYCOPY_H_
#define UU_OPERATIONS_EMPTYCOPY_H_

#include <memory>
#include <string>
#include "networks/Network.hpp"
#include "networks/ProbabilisticNetwork.hpp"

namespace uu {
namespace net {

/**
 * Returns an empty graph with the same features as G (directed/undirected, allowing loops, ...).
 */
template<typename G>
std::unique_ptr<G>
empty_copy(
    const G* g,
    const std::string& name = ""
);

}
}

#include "empty_copy.ipp"

#endif
