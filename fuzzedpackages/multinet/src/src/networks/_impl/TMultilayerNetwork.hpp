/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_MNET_DATASTRUCTURE_GRAPHS_MULTILAYERNETWORK_H_
#define UU_MNET_DATASTRUCTURE_GRAPHS_MULTILAYERNETWORK_H_

#include <memory>
#include <string>
#include <unordered_set>
#include "core/exceptions/ElementNotFoundException.hpp"
#include "core/datastructures/observers/ObserverStore.hpp"
#include "objects/Edge.hpp"
#include "objects/Vertex.hpp"
#include "networks/_impl/TMultilayerNetworkType.hpp"

namespace uu {
namespace net {

/**
 * TMultilayerNetwork is a generic class that can be instantiated into several specific types of
 * network whose vertices/edges are organized into layers.
 *
 * A TMultilayerNetwork is instantiated by specifying a VertexStore (V), containing all the vertices
 * in the network, a LayerStore (L) handling the individual layers (each corresponding to a graphs
 * ever the vertices in V), and an EdgeStore (E) for inter-layer edges.
 */
template <typename V, typename L, typename E>
class TMultilayerNetwork
    :
    public core::ObserverStore
{

  public:

    /**
     * Creates an empty network.
     */
    TMultilayerNetwork(
        const std::string& name,
        TMultilayerNetworkType t,
        std::unique_ptr<V> v,
        std::unique_ptr<L> l,
        std::unique_ptr<E> e
    );

    /**
     * Returns a pointer to the network's vertex store.
     */
    V*
    vertices(
    );

    /**
     * Returns a pointer to the network's (const) vertex store.
     */
    const V*
    vertices(
    ) const;


    /**
     * Returns a pointer to the network's edge store.
     */
    L*
    layers(
    );


    /**
     * Returns a pointer to the network's (const) edge store.
     */
    const L*
    layers(
    ) const;

    bool
    is_ordered(
    ) const;

    /**
     * Returns a string providing a summary of the graph structure.
     */
    std::string
    summary(
    ) const;

  public:

    const std::string name;

    /**
     * Returns a pointer to the container of edge stores for each pair of layers.
     */
    E*
    interlayer_edges(
    );


    /**
     * Returns a pointer to the container of edge stores for each pair of layers.
     */
    const E*
    interlayer_edges(
    ) const;

  protected:

    /** Internal vertex store. */
    std::unique_ptr<V> vertices_;

    /** Internal layer store. */
    std::unique_ptr<L> layers_;

    /** Internal edge store. */
    std::unique_ptr<E> edges_;

    TMultilayerNetworkType type_;

    /** ... */
    std::unordered_set<std::unique_ptr<core::GenericObserver>> obs_;

};



template <typename V, typename L, typename E>
TMultilayerNetwork<V,L,E>::
TMultilayerNetwork(
    const std::string& name,
    TMultilayerNetworkType t,
    std::unique_ptr<V> v,
    std::unique_ptr<L> l,
    std::unique_ptr<E> e
) : name(name)
{
    vertices_ = std::move(v);
    layers_ = std::move(l);
    edges_ = std::move(e);
    type_ = t;
}


template <typename V, typename L, typename E>
V*
TMultilayerNetwork<V,L,E>::
vertices(
)
{
    return vertices_.get();
}


template <typename V, typename L, typename E>
const V*
TMultilayerNetwork<V,L,E>::
vertices(
) const
{
    return vertices_.get();
}


template <typename V, typename L, typename E>
E*
TMultilayerNetwork<V,L,E>::
interlayer_edges(
)
{
    return edges_.get();
}


template <typename V, typename L, typename E>
const E*
TMultilayerNetwork<V,L,E>::
interlayer_edges(
) const
{
    return edges_.get();
}


template <typename V, typename L, typename E>
L*
TMultilayerNetwork<V,L,E>::
layers(
)
{
    return layers_.get();
}


template <typename V, typename L, typename E>
const L*
TMultilayerNetwork<V,L,E>::
layers(
) const
{
    return layers_.get();
}


template <typename V, typename L, typename E>
bool
TMultilayerNetwork<V,L,E>::
is_ordered(
) const
{
    return type_.is_ordered;
}

template <typename V, typename L, typename E>
std::string
TMultilayerNetwork<V,L,E>::
summary(
) const
{
    std::string summary =
        "TMultilayerNetwork (" +
        vertices_->summary() + ", " +
        layers_->summary() + ", " +
        edges_->summary() + ")";
    return summary;
}


}
}

#endif
