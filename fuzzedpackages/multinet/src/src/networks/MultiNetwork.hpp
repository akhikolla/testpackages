#ifndef UU_NETWORKS_MULTINETWORK_H_
#define UU_NETWORKS_MULTINETWORK_H_

#include <memory>
#include <string>
#include "networks/_impl/Graph.hpp"
#include "networks/_impl/stores/AttrVertexStore.hpp"
#include "networks/_impl/stores/AttrMultiEdgeStore.hpp"

namespace uu {
namespace net {

/**
 * A MultiNetwork is an attributed graph allowing multiple edges between each pair of vertices.
 *
 * Vertex and edge attributes are local to the network, that is, the same vertex inside another
 * network will have different attributes.
 *
 * Depending on its parameters, a MultiNetwork can allow or disallow loops (default: allow) and
 * be directed or undirected (default: undirected). That is, a MultiNetwork by default corresponds
 * to a mathematical multigraph.
 */
class MultiNetwork
{

  public:

    const std::string name;

    /**
     * Creates a MultiNetwork with directed or undirected multiedges and with or without loops.
     */
    MultiNetwork(
        const std::string& name,
        EdgeDir dir = EdgeDir::UNDIRECTED,
        bool allow_loops = true
    );

    /**
     * Returns a pointer to the network's vertices.
     */
    AttrVertexStore*
    vertices(
    );

    /**
     * Returns a pointer to the network's (const) vertices.
     */
    const AttrVertexStore*
    vertices(
    ) const;


    /**
     * Returns a pointer to the network's edges.
     */
    AttrMultiEdgeStore*
    edges(
    );


    /**
     * Returns a pointer to the network's (const) edges.
     */
    const AttrMultiEdgeStore*
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
     * Checks if the network has temporal information on its edges.
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
     * Always returns true.
     */
    virtual
    bool
    allows_multi_edges(
    ) const;


  private:

    std::unique_ptr<Graph<AttrVertexStore, AttrMultiEdgeStore>> data_;

};

}
}

#endif
