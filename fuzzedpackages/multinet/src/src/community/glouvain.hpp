#ifndef UU_COMMUNITY_GLOUVAIN_H_
#define UU_COMMUNITY_GLOUVAIN_H_


#include "core/utils/pretty_printing.hpp"
#include "community/CommunityStructure.hpp"
#include "community/VertexLayerCommunity.hpp"
#include "community/_impl/cutils.hpp"
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include <list>
#include <limits>

namespace uu {
namespace net {

template <typename M, typename G>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const G>>>
generalized_louvain(
    const M* mnet,
    double gamma,
    double omega,
    size_t limit
);


}
}

#include "glouvain.ipp"

#endif
