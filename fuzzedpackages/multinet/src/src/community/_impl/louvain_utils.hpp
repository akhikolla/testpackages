#ifndef UU_NET_COMMUNITY_LOUVAINUTILS_H_
#define UU_NET_COMMUNITY_LOUVAINUTILS_H_

#include <chrono>
#include <unordered_map>
#include <vector>
#include <memory>
#include "community/_impl/MetaNetwork.hpp"
#include "community/CommunityStructure.hpp"
#include "community/Community.hpp"
#include "objects/EdgeMode.hpp"
#include "objects/Vertex.hpp"
#include "measures/size.hpp"
#include "measures/strength.hpp"

namespace uu {
namespace net {


std::unique_ptr<MetaNetwork>
convert(
    const Network* g
);


std::unique_ptr<MetaNetwork>
aggregate(
    const WeightedNetwork* g,
    std::unordered_map<const Vertex*, size_t> community
);


void
expand(
    const std::vector<std::unique_ptr<MetaNetwork>>& levels,
    size_t i,
    const Vertex* v,
    Community<const Vertex*>* com
);

std::unique_ptr<CommunityStructure<Community<const Vertex*>>>
communities(
    const std::vector<std::unique_ptr<MetaNetwork>>& levels
);

std::unique_ptr<MetaNetwork>
pass(
    const WeightedNetwork* g
);


}
}

#endif
