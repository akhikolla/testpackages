#ifndef UU_NETWORKS_IMPL_STORES_MLSIMPLEEDGESTORE_H_
#define UU_NETWORKS_IMPL_STORES_MLSIMPLEEDGESTORE_H_

#include "networks/_impl/stores/AttributedDynamicInterlayerSimpleEdgeStore.hpp"
#include "objects/Vertex.hpp"
#include "networks/Network.hpp"
#include "networks/_impl/stores/Attributes.hpp"
#include "networks/_impl/stores/UserDefinedAttrs.hpp"


namespace uu {
namespace net {

typedef  AttributedDynamicInterlayerSimpleEdgeStore<Vertex,Network,Attributes<MLEdge<Vertex,Network>, UserDefinedAttrs<MLEdge<Vertex,Network>>>> MLSimpleEdgeStore;

}
}

#endif
