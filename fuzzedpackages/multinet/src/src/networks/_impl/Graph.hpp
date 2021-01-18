/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_DATASTRUCTURE_GRAPHS_GRAPH_H_
#define UU_NET_DATASTRUCTURE_GRAPHS_GRAPH_H_

#include <memory>
#include <string>
#include "core/exceptions/ElementNotFoundException.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include "core/datastructures/observers/ObserverStore.hpp"
#include "networks/_impl/observers/AdjVertexCheckObserver.hpp"
#include "networks/_impl/observers/PropagateObserver.hpp"
#include "networks/_impl/GraphType.hpp"

namespace uu {
namespace net {

/**
 * Graph is a generic class that can be instantiated into several specific types of network.
 *
 * A Graph is instantiated specifying a VertexStore (of type V) and an EdgeStore (of type E).
 * GT is an object of type GraphType, that provides a number of predicates in the form is_xxx()
 * that can be used by functions to check whether the graph implements specific features such as
 * being weighted, allowing user-defined attributes, etc.
 */
template<typename V, typename E>
class Graph
    :
    public core::ObserverStore
{

  public:

    /**
     * Creates an empty graph.
     *
     * This is only used internally. Use:
     * - Graph::create() to create a new graph,
     * - an IO function available in the library to read a graph from file, or
     * - a graph generation function to generate a synthetic graph.
     *
     */
    Graph(
        const std::string& name,
        GraphType t,
        std::unique_ptr<V> v,
        std::unique_ptr<E> e
    );

    /**
     * Returns a pointer to the graph's vertex store.
     */
    V*
    vertices(
    );

    /**
     * Returns a pointer to the graph's (const) vertex store.
     */
    const V*
    vertices(
    ) const;


    /**
     * Returns a pointer to the graph's edge store.
     */
    E*
    edges(
    );


    /**
     * Returns a pointer to the graph's (const) edge store.
     */
    const E*
    edges(
    ) const;


    /**
     * Checks if the edges in this graph are directed.
     */
    bool
    is_directed(
    ) const;


    /**
     * Checks if the graph is weighted.
     */
    bool
    is_weighted(
    ) const;

    /**
     * Checks if the graph is probabilistic.
     */
    bool
    is_probabilistic(
    ) const;


    /**
     * Checks if the graph has temporal information on its edges.
     */
    bool
    is_temporal(
    ) const;


    /**
     * Checks if the graph allows users to define their own generic attributes.
     */
    bool
    is_attributed(
    ) const;


    /**
     * Checks if the graph allows multi-edges. If false, only simple edges are allowed.
     */
    bool
    allows_multi_edges(
    ) const;


    /**
     * Checks if the graph allows loops.
     */
    bool
    allows_loops(
    ) const;


    /**
     * Returns a string providing a summary of the graph structure.
     */
    virtual
    std::string
    summary(
    ) const;

  public:

    const std::string name;

  protected:

    /** Graph type. */
    GraphType type_;

    /** Internal vertex store. */
    std::unique_ptr<V> vertices_;

    /** Internal edge store. */
    std::unique_ptr<E> edges_;

};



template<typename V, typename E>
Graph<V,E>::
Graph(
    const std::string& name,
    GraphType t,
    std::unique_ptr<V> v,
    std::unique_ptr<E> e
) : name(name), type_(t)
{
    vertices_ = std::move(v);
    edges_ = std::move(e);

    if (edges_->is_directed() != t.is_directed)
    {
        throw core::WrongParameterException("incompatible graph type directionality and edge store directionality");
    }


    // register an observer to propagate the removal of vertices to the edge store
    auto obs1 = std::make_unique<PropagateObserver<E, const Vertex>>(edges());
    vertices()->attach(obs1.get());
    register_observer(std::move(obs1));

    // register an observer to check that the end vertices of a newly inserted graph exist
    auto obs2 = std::make_unique<AdjVertexCheckObserver<V>>(vertices());
    edges()->attach(obs2.get());
    register_observer(std::move(obs2));

}

template<typename V, typename E>
V*
Graph<V,E>::
vertices(
)
{
    return vertices_.get();
}


template<typename V, typename E>
const V*
Graph<V,E>::
vertices(
) const
{
    return vertices_.get();
}


template<typename V, typename E>
E*
Graph<V,E>::
edges(
)
{
    return edges_.get();
}


template<typename V, typename E>
const E*
Graph<V,E>::
edges(
) const
{
    return edges_.get();
}


template<typename V, typename E>
bool
Graph<V,E>::
is_directed(
) const
{
    return type_.is_directed;
}


template<typename V, typename E>
bool
Graph<V,E>::
is_weighted() const
{
    return type_.is_weighted;
}


template<typename V, typename E>
bool
Graph<V,E>::
is_probabilistic() const
{
    return type_.is_probabilistic;
}


template<typename V, typename E>
bool
Graph<V,E>::
is_temporal() const
{
    return type_.is_temporal;
}

template<typename V, typename E>
bool
Graph<V,E>::
is_attributed() const
{
    return type_.is_attributed;
}

template<typename V, typename E>
bool
Graph<V,E>::
allows_multi_edges() const
{
    return type_.allows_multi_edges;
}

template<typename V, typename E>
bool
Graph<V,E>::
allows_loops() const
{
    return type_.allows_loops;
}


template<typename V, typename E>
std::string
Graph<V,E>::
summary(
) const
{
    std::string summary =
        "Graph (" +
        vertices_->summary() + ", " +
        edges_->summary() + ")";
    return summary;
}


}
}

#endif
