#ifndef UU_COMMUNITY_GLOUVAIN2_H_
#define UU_COMMUNITY_GLOUVAIN2_H_


#include "community/CommunityStructure.hpp"
#include "community/VertexLayerCommunity.hpp"
#include "objects/Vertex.hpp"
#include "community/_impl/glouvain2_utils.hpp"
#include <memory>


namespace uu {
namespace net {

template <typename M>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const typename M::layer_type>>>
glouvain2(
    const M* g,
    double omega
);

}
}

#include "glouvain2.ipp"

#endif
