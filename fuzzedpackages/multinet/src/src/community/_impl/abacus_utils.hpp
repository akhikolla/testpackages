#ifndef UU_COMMUNITY_IMPL_ABACUSUTILS_H_
#define UU_COMMUNITY_IMPL_ABACUSUTILS_H_

#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <random>
#include "community/CommunityStructure.hpp"
#include "community/VertexLayerCommunity.hpp"
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include <vector>
#include "community/Community.hpp"
#include "objects/EdgeMode.hpp"
#include "objects/Vertex.hpp"

extern "C" {
#include <eclat.h>
}

namespace uu {
namespace net {

/*  */
template <typename M, typename L>
void
read_layers(
    const M* mnet,
    PillarCommunity<L>* com,
    FILE* file
);

/*  */
template <typename M, typename L>
int
read_actors(
    const M* mnet,
    PillarCommunity<L>* com,
    FILE* tidfile
);

/*  */
template <typename M, typename L>
std::unique_ptr<CommunityStructure<PillarCommunity<L>>>
read_eclat_communities(
    const M* mnet,
    FILE* file,
    FILE* tidfile
);

/*  */
template <typename M, typename L>
std::unique_ptr<CommunityStructure<PillarCommunity<L>>>
eclat_merge(
    const M* mnet,
    const std::unordered_map<const L*, CommunityStructure<Community<const Vertex*>>*>& single_layer_communities,
    int min_actors,
    int min_layers
);


}
}

#include "abacus_utils.ipp"

#endif
