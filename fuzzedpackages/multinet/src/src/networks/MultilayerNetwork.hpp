#ifndef UU_NETWORKS_MULTILAYERNETWORK_H_
#define UU_NETWORKS_MULTILAYERNETWORK_H_

/**
 *
 */

#include <memory>
#include <string>
#include "networks/_impl/TMultilayerNetwork.hpp"
#include "networks/_impl/stores/AttrVertexStore.hpp"
#include "networks/_impl/stores/MLLayerStore.hpp"
#include "networks/_impl/stores/MLSimpleEdgeStore.hpp"

#include "objects/Vertex.hpp"
#include "networks/Network.hpp"

namespace uu {
namespace net {

/**
 * A MultilayerNetwork is a set of Networks (called layers) allowing edges across layers.
 *
 * Vertices are called actors, and the same actor can be present in multiple layers.
 * A multilayer vertex (MLVertex) is a pair (actor,layer).
 * A multilayer edge (MLEdge) is an edge connecting multilayer vertices.
 */
class
    MultilayerNetwork
{

  public:

    const std::string name;

    typedef Network layer_type;
    typedef Vertex vertex_type;

    /**
     * Creates a Network with no layers.
     */
    MultilayerNetwork(
        const std::string& name
    );


    /**
     * Returns a pointer to the network's actors.
     */
    AttrVertexStore*
    actors(
    );


    /**
     * Returns a pointer to the network's actors.
     */
    const AttrVertexStore*
    actors(
    ) const;


    /**
     * Returns a pointer to the network's layers.
     */
    MLLayerStore*
    layers(
    );


    /**
     * Returns a pointer to the network's layers.
     */
    const MLLayerStore*
    layers(
    ) const;


    /**
     * Returns a pointer to the network's interlayer edges.
     */
    MLSimpleEdgeStore*
    interlayer_edges(
    );


    /**
     * Returns a pointer to the network's interlayer edges.
     */
    const MLSimpleEdgeStore*
    interlayer_edges(
    ) const;


    /**
     * Checks if the network allows interlayer edges.
     * Always returns false.
     */
    bool
    is_ordered(
    ) const;


    /**
     * Checks if the network allows interlayer edges.
     * Always returns true.
     */
    bool
    allows_interlayer_edges(
    ) const;


    std::string
    summary(
    ) const;


  private:

    std::unique_ptr<TMultilayerNetwork<AttrVertexStore, MLLayerStore, MLSimpleEdgeStore>> data_;

};


}
}

#endif
