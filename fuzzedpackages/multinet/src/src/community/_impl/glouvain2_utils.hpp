#ifndef UU_NET_COMMUNITY_GLOUVAINUTILS_H_
#define UU_NET_COMMUNITY_GLOUVAINUTILS_H_

#include <chrono>
#include <unordered_map>
#include <vector>
#include <memory>
#include "community/_impl/GMetaNetwork.hpp"
#include "community/CommunityStructure.hpp"
#include "community/Community.hpp"
#include "objects/EdgeMode.hpp"
#include "objects/Vertex.hpp"
#include "measures/size.hpp"
#include "measures/strength.hpp"
#include "networks/MultilayerNetwork.hpp"
#include "networks/OrderedMultiplexNetwork.hpp"

namespace uu {
namespace net {


std::tuple<std::unique_ptr<GMetaNetwork>, std::map<const Vertex*, std::pair<const Vertex*, const Network*>>, std::vector<std::unique_ptr<const Vertex>>>
        convert(
            const MultilayerNetwork* g,
            double omega
        );

std::tuple<std::unique_ptr<GMetaNetwork>, std::map<const Vertex*, std::pair<const Vertex*, const Network*>>, std::vector<std::unique_ptr<const Vertex>>>
        convert(
            const OrderedMultiplexNetwork* g,
            double omega
        );

void
expand(
    const std::vector<std::unique_ptr<GMetaNetwork>>& levels,
    size_t i,
    const Vertex* v,
    Community<const Vertex*>* com
);

std::unique_ptr<CommunityStructure<Community<const Vertex*>>>
communities(
    const std::vector<std::unique_ptr<GMetaNetwork>>& levels
);

std::unique_ptr<GMetaNetwork>
aggregate(
    const GMetaNetwork* meta,
    std::unordered_map<const Vertex*, size_t> community
);

std::unique_ptr<GMetaNetwork>
pass(
    const GMetaNetwork* meta
);


}
}

#endif
