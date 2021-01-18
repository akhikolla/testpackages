#ifndef UU_MNET_DATASTRUCTURE_GRAPHS_ATTRIBUTEDORDEREDHOMOGENEOUSMULTILAYERNETWORK_H_
#define UU_MNET_DATASTRUCTURE_GRAPHS_ATTRIBUTEDORDEREDHOMOGENEOUSMULTILAYERNETWORK_H_

#include <memory>
#include <string>
#include "networks/_impl/TMultilayerNetwork.hpp"
#include "networks/_impl/stores/AttrVertexStore.hpp"
#include "networks/_impl/stores/MLOrderedLayerStore.hpp"
#include "networks/_impl/stores/MLSimpleEdgeStore.hpp"

#include "objects/Vertex.hpp"
#include "networks/Network.hpp"

namespace uu {
namespace net {

/**
 * An OrderedMultiplexNetwork is a set of ordered Networks (called layers).
 *
 * Vertices are called actors, and the same actor can be present in multiple layers.
 * A multilayer vertex (MLVertex) is a pair (actor,layer).
 */
class
    OrderedMultiplexNetwork
{

  public:

    const std::string name;

    typedef Network layer_type;
    typedef Vertex vertex_type;

    OrderedMultiplexNetwork(
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
    MLOrderedLayerStore*
    layers(
    );


    /**
     * Returns a pointer to the network's layers.
     */
    const MLOrderedLayerStore*
    layers(
    ) const;


    /**
     * Checks if the network allows interlayer edges.
     * Always returns false.
     */
    bool
    is_ordered(
    ) const
    {
        return false;
    }


    /**
     * Checks if the network allows interlayer edges.
     * Always returns true.
     */
    bool
    allows_interlayer_edges(
    ) const
    {
        return false;
    }


    std::string
    summary(
    ) const;


  private:

    std::unique_ptr<TMultilayerNetwork<AttrVertexStore, MLOrderedLayerStore, MLSimpleEdgeStore>> data_;

};

}
}

#endif
