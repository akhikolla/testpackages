#ifndef UU_NETWORKS_NETWORK_H_
#define UU_NETWORKS_NETWORK_H_

#include <memory>
#include <string>
#include "networks/_impl/Graph.hpp"
#include "networks/_impl/stores/AttrVertexStore.hpp"
#include "networks/_impl/stores/AttrSimpleEdgeStore.hpp"

namespace uu {
namespace net {

/**
 * A Network is an attributed graph with at most one edge between each pair of vertices.
 *
 * Vertex and edge attributes are local to the network, that is, the same vertex inside another
 * network will have different attributes.
 * Depending on its parameters, a Network can allow or disallow loops (default: disallow) and
 * be directed or undirected (default: undirected). That is, a Network by default corresponds to
 * a mathematical simple graph.
 */
class Network
{

  public:

    const std::string name;

    /**
     * Creates a Network with directed or undirected simple edges and allowing or not loops.
     */
    Network(
        const std::string& name,
        EdgeDir dir = EdgeDir::UNDIRECTED,
        bool allow_loops = false
    );


    /**
     * Returns a pointer to the network's vertices.
     */
    AttrVertexStore*
    vertices(
    );


    /**
     * Returns a pointer to the network's vertices.
     */
    const AttrVertexStore*
    vertices(
    ) const;


    /**
     * Returns a pointer to the network's edges.
     */
    AttrSimpleEdgeStore*
    edges(
    );


    /**
     * Returns a pointer to the network's edges.
     */
    const AttrSimpleEdgeStore*
    edges(
    ) const;


    /**
     * Checks if the edges in this network are directed.
     */
    virtual
    bool
    is_directed(
    ) const;


    /**
     * Checks if the network allows loops.
     */
    virtual
    bool
    allows_loops(
    ) const;


    /**
     * Checks if the network is weighted.
     * Always returns false.
     */
    virtual
    bool
    is_weighted(
    ) const;


    /**
     * Checks if the network is probabilistic.
     * Always returns false.
     */
    virtual
    bool
    is_probabilistic(
    ) const;


    /**
     * Checks if the network is temporal.
     * Always returns false.
     */
    virtual
    bool
    is_temporal(
    ) const;


    /**
     * Checks if the network allows users to define their own generic attributes.
     * Always returns true.
     */
    virtual
    bool
    is_attributed(
    ) const;


    /**
     * Checks if the network allows multi-edges.
     * Always returns false.
     */
    virtual
    bool
    allows_multi_edges(
    ) const;


  private:

    std::unique_ptr<Graph<AttrVertexStore, AttrSimpleEdgeStore>> data_;

};

}
}

#endif
