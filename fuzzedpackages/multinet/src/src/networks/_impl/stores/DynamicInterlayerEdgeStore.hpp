/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_DATASTRUCTURES_STORES_DYNAMICINTERLAYEREDGESTORE_H_
#define UU_NET_DATASTRUCTURES_STORES_DYNAMICINTERLAYEREDGESTORE_H_

#include <memory>
#include <unordered_map>
#include <unordered_map>
#include "core/datastructures/containers/SharedPtrSortedRandomSet.hpp"
#include "core/datastructures/observers/Subject.hpp"
#include "objects/MLEdge.hpp"
#include "objects/EdgeDir.hpp"
#include "objects/EdgeMode.hpp"
#include "networks/_impl/containers/GenericObjectList.hpp"
#include "core/exceptions/OperationNotSupportedException.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace net {

template <typename V, typename L>
class DynamicInterlayerEdgeStore:
    public core::SharedPtrSortedRandomSet<const MLEdge<V,L>>,
            public core::Subject<const MLEdge<V,L>>
{

  private:

    typedef core::SharedPtrSortedRandomSet<const MLEdge<V,L>> super;

  public:

    /**
     * Constructor.
     */

    //DynamicInterlayerEdgeStore(
    //);

  public:

    using super::size;
    using super::add;
    using super::erase;

    virtual
    const MLEdge<V,L>*
    add(
        std::shared_ptr<const MLEdge<V,L>> e
    ) override;

    /**
     * Adds a new edge.
     * Multiple edges between the same pair of vertices are not allowed.
     * @param vertex1 a pointer to the "from" vertex if directed, or to one end of
     * the edge if undirected.
     * @param vertex2 a pointer to the "to" vertex if directed, or one end of the
     * edge if undirected.
     * @return a pointer to the new edge, or nullptr if the edge already exists.
     **/
    virtual
    const MLEdge<V,L> *
    add(
        const V* vertex1,
        const L* layer1,
        const V* vertex2,
        const L* layer2
    );

    virtual
    bool
    erase(
        const MLEdge<V,L>* e
    ) override = 0;


    virtual
    GenericObjectList<MLEdge<V,L>>*
                                get(
                                    const L* layer1,
                                    const L* layer2
                                ) const;

    /**
     * @brief Returns the nodes with an edge from/to the input vertex.
     * @param node pointer to the node.
     * @param mode IN, OUT or INOUT.
     * @return the list of neighbors.
     **/
    const
    GenericObjectList<V>*
    neighbors(
        const L* layer1,
        const L* layer2,
        const V* vertex,
        EdgeMode mode
    ) const;

    /**
     * @brief Returns the nodes with an edge from/to the input vertex.
     * @param node pointer to the node.
     * @param mode IN, OUT or INOUT.
     * @return the list of neighbors.
     **/
    const
    GenericObjectList<MLEdge<V,L>>*
                                incident(
                                    const L* layer1,
                                    const L* layer2,
                                    const V* vertex,
                                    EdgeMode mode
                                ) const;


    bool
    is_directed(
        const L* layer1,
        const L* layer2
    ) const;


    void
    set_directed(
        const L* layer1,
        const L* layer2,
        bool directed
    );


    virtual
    void
    add(
        const L* layer
    );

    virtual
    void
    erase(
        const L* layer
    );


    virtual
    void
    erase(
        const L* layer,
        const V* vertex
    ) = 0;

  protected:


    /** Edges */
    std::unordered_map<const L*, std::unordered_map<const L*, std::unique_ptr<GenericObjectList<MLEdge<V,L>>>>> edges_;

    /** Edge directionality */
    std::unordered_map<const L*, std::unordered_map<const L*, EdgeDir>> edge_directionality;

    // Indexes to sets of objects (Set IDX):
    std::unordered_map<const L*, std::unordered_map<const L*, std::unordered_map<const V*, std::unique_ptr<GenericObjectList<V>>>>> sidx_neighbors_out;
    std::unordered_map<const L*, std::unordered_map<const L*, std::unordered_map<const V*, std::unique_ptr<GenericObjectList<V>>>>> sidx_neighbors_in;
    std::unordered_map<const L*, std::unordered_map<const L*, std::unordered_map<const V*, std::unique_ptr<GenericObjectList<V>>>>> sidx_neighbors_all;

    std::unordered_map<const L*, std::unordered_map<const L*, std::unordered_map<const V*, std::unique_ptr<GenericObjectList<MLEdge<V,L>>>>>> sidx_incident_out;
    std::unordered_map<const L*, std::unordered_map<const L*, std::unordered_map<const V*, std::unique_ptr<GenericObjectList<MLEdge<V,L>>>>>> sidx_incident_in;
    std::unordered_map<const L*, std::unordered_map<const L*, std::unordered_map<const V*, std::unique_ptr<GenericObjectList<MLEdge<V,L>>>>>> sidx_incident_all;
};

template <typename V, typename L>
const MLEdge<V,L> *
DynamicInterlayerEdgeStore<V,L>::
add(
    const V* vertex1,
    const L* layer1,
    const V* vertex2,
    const L* layer2
)
{
    core::assert_not_null(vertex1, "add", "vertex1");
    core::assert_not_null(layer1, "add", "layer1");
    core::assert_not_null(vertex2, "add", "vertex2");
    core::assert_not_null(layer2, "add", "layer2");

    auto dir = edge_directionality.at(layer1).at(layer2);
    auto edge = MLEdge<V,L>::create(vertex1, layer1, vertex2, layer2, dir);
    return add(edge);
}


template <typename V, typename L>
const MLEdge<V,L>*
DynamicInterlayerEdgeStore<V,L>::
add(
    std::shared_ptr<const MLEdge<V,L>> e
)
{
    core::assert_not_null(e.get(), "add", "e");

    if (e->dir != edge_directionality.at(e->l1).at(e->l2))
    {
        throw core::OperationNotSupportedException("wrong edge directionality");
    }

    const MLEdge<V,L>* new_edge = super::add(e);

    if (!new_edge) // edge already existing
    {
        return nullptr;
    }

    edges_.at(e->l1).at(e->l2)->add(new_edge);
    edges_.at(e->l2).at(e->l1)->add(new_edge);

    if (sidx_neighbors_out[e->l1][e->l2].count(e->v1)==0)
    {
        sidx_neighbors_out[e->l1][e->l2][e->v1] = std::make_unique<GenericObjectList<V>>();
        sidx_incident_out[e->l1][e->l2][e->v1] = std::make_unique<GenericObjectList<MLEdge<V,L>>>();
    }

    sidx_neighbors_out[e->l1][e->l2][e->v1]->add(e->v2);
    sidx_incident_out[e->l1][e->l2][e->v1]->add(new_edge);


    if (sidx_neighbors_in[e->l2][e->l1].count(e->v2)==0)
    {
        sidx_neighbors_in[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<V>>();
        sidx_incident_in[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<MLEdge<V,L>>>();
    }

    sidx_neighbors_in[e->l2][e->l1][e->v2]->add(e->v1);
    sidx_incident_in[e->l2][e->l1][e->v2]->add(new_edge);


    if (sidx_neighbors_all[e->l1][e->l2].count(e->v1)==0)
    {
        sidx_neighbors_all[e->l1][e->l2][e->v1] = std::make_unique<GenericObjectList<V>>();
        sidx_incident_all[e->l1][e->l2][e->v1] = std::make_unique<GenericObjectList<MLEdge<V,L>>>();
    }

    sidx_neighbors_all[e->l1][e->l2][e->v1]->add(e->v2);
    sidx_incident_all[e->l1][e->l2][e->v1]->add(new_edge);

    if (sidx_neighbors_all[e->l2][e->l1].count(e->v2)==0)
    {
        sidx_neighbors_all[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<V>>();
        sidx_incident_all[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<MLEdge<V,L>>>();
    }

    sidx_neighbors_all[e->l2][e->l1][e->v2]->add(e->v1);
    sidx_incident_all[e->l2][e->l1][e->v2]->add(new_edge);


    if (e->dir == EdgeDir::UNDIRECTED)
    {

        if (sidx_neighbors_out[e->l2][e->l1].count(e->v2)==0)
        {
            sidx_neighbors_out[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<V>>();
            sidx_incident_out[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<MLEdge<V,L>>>();
        }

        sidx_neighbors_out[e->l2][e->l1][e->v2]->add(e->v1);
        sidx_incident_out[e->l2][e->l1][e->v2]->add(new_edge);

        if (sidx_neighbors_in[e->l1][e->l2].count(e->v1)==0)
        {
            sidx_neighbors_in[e->l1][e->l2][e->v1] = std::make_unique<GenericObjectList<V>>();
            sidx_incident_in[e->l1][e->l2][e->v1] = std::make_unique<GenericObjectList<MLEdge<V,L>>>();
        }

        sidx_neighbors_in[e->l1][e->l2][e->v1]->add(e->v2);
        sidx_incident_in[e->l1][e->l2][e->v1]->add(new_edge);

        /*
        if (sidx_neighbors_all[e->l2][e->l1].count(e->v2)==0)
        {
            sidx_neighbors_all[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<V>>();
            sidx_incident_all[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<MLEdge<V,L>>>();
        }

        sidx_neighbors_all[e->l2][e->l1][e->v2]->add(e->v1);
        sidx_incident_all[e->l2][e->l1][e->v2]->add(new_edge);
         */
    }

    return new_edge;
}

template <typename V, typename L>
GenericObjectList<MLEdge<V,L>>*
                            DynamicInterlayerEdgeStore<V,L>::
                            get(
                                const L* layer1,
                                const L* layer2
                            ) const
{
    core::assert_not_null(layer1, "neighbors", "layer1");
    core::assert_not_null(layer2, "neighbors", "layer2");
    return edges_.at(layer1).at(layer2).get();
}

template <typename V, typename L>
const
GenericObjectList<V>*
DynamicInterlayerEdgeStore<V,L>::
neighbors(
    const L* layer1,
    const L* layer2,
    const V* vertex,
    EdgeMode mode
) const
{

    core::assert_not_null(layer1, "neighbors", "layer1");
    core::assert_not_null(layer2, "neighbors", "layer2");
    core::assert_not_null(vertex, "neighbors", "vertex");

    if (mode==EdgeMode::IN)
    {
        if (sidx_neighbors_in.at(layer1).at(layer2).count(vertex)==0)
        {
            return GenericObjectList<V>::empty.get();
        }

        return sidx_neighbors_in.at(layer1).at(layer2).at(vertex).get();
    }

    else if (mode==EdgeMode::OUT)
    {
        if (sidx_neighbors_out.at(layer1).at(layer2).count(vertex)==0)
        {
            return GenericObjectList<V>::empty.get();
        }

        return sidx_neighbors_out.at(layer1).at(layer2).at(vertex).get();
    }

    else if (mode==EdgeMode::INOUT)
    {
        if (sidx_neighbors_all.at(layer1).at(layer2).count(vertex)==0)
        {
            return GenericObjectList<V>::empty.get();
        }

        return sidx_neighbors_all.at(layer1).at(layer2).at(vertex).get();
    }

    else
    {
        throw core::WrongParameterException("neighborhood mode");
    }
}


template <typename V, typename L>
const
GenericObjectList<MLEdge<V,L>>*
                            DynamicInterlayerEdgeStore<V,L>::
                            incident(
                                const L* layer1,
                                const L* layer2,
                                const V* vertex,
                                EdgeMode mode
                            ) const
{

    core::assert_not_null(layer1, "neighbors", "layer1");
    core::assert_not_null(layer2, "neighbors", "layer2");
    core::assert_not_null(vertex, "incident", "vertex");

    if (mode==EdgeMode::IN)
    {
        if (sidx_incident_in.at(layer1).at(layer2).count(vertex)==0)
        {
            return GenericObjectList<MLEdge<V,L>>::empty.get();
        }

        return sidx_incident_in.at(layer1).at(layer2).at(vertex).get();
    }

    else if (mode==EdgeMode::OUT)
    {
        if (sidx_incident_out.at(layer1).at(layer2).count(vertex)==0)
        {
            return GenericObjectList<MLEdge<V,L>>::empty.get();
        }

        return sidx_incident_out.at(layer1).at(layer2).at(vertex).get();
    }

    else if (mode==EdgeMode::INOUT)
    {
        if (sidx_incident_all.at(layer1).at(layer2).count(vertex)==0)
        {
            return GenericObjectList<MLEdge<V,L>>::empty.get();
        }

        return sidx_incident_all.at(layer1).at(layer2).at(vertex).get();
    }

    else
    {
        throw core::WrongParameterException("neighborhood mode");
    }
}


template <typename V, typename L>
bool
DynamicInterlayerEdgeStore<V,L>::
is_directed(
    const L* layer1,
    const L* layer2
) const
{
    core::assert_not_null(layer1, "is_directed", "layer1");
    core::assert_not_null(layer2, "is_directed", "layer2");
    return edge_directionality.at(layer1).at(layer2) == EdgeDir::DIRECTED?true:false;
}


template <typename V, typename L>
void
DynamicInterlayerEdgeStore<V,L>::
set_directed(
    const L* layer1,
    const L* layer2,
    bool directed
)
{
    core::assert_not_null(layer1, "set_directed", "layer1");
    core::assert_not_null(layer2, "set_directed", "layer2");

    if (edges_.at(layer1).at(layer2)->size() > 0)
    {
        throw core::OperationNotSupportedException("cannot change directionality after edges have been inserted");
    }

    edge_directionality.at(layer1).at(layer2) = directed?EdgeDir::DIRECTED:EdgeDir::UNDIRECTED;
    edge_directionality.at(layer2).at(layer1) = directed?EdgeDir::DIRECTED:EdgeDir::UNDIRECTED;
}



template <typename V, typename L>
void
DynamicInterlayerEdgeStore<V,L>::
add(
    const L* layer
)
{
    core::assert_not_null(layer, "add", "layer");

    edges_[layer];

    std::vector<const L*> layers;

    for (auto&& p: edges_)
    {
        layers.push_back(p.first);
    }

    for (auto l: layers)
    {
        edges_[l][layer] = std::make_unique<GenericObjectList<MLEdge<V,L>>>();
        edge_directionality[l][layer];
        sidx_neighbors_out[l][layer];
        sidx_neighbors_in[l][layer];
        sidx_neighbors_all[l][layer];
        sidx_incident_out[l][layer];
        sidx_incident_in[l][layer];
        sidx_incident_all[l][layer];

        edges_[layer][l] = std::make_unique<GenericObjectList<MLEdge<V,L>>>();
        edge_directionality[layer][l];
        sidx_neighbors_out[layer][l];
        sidx_neighbors_in[layer][l];
        sidx_neighbors_all[layer][l];
        sidx_incident_out[layer][l];
        sidx_incident_in[layer][l];
        sidx_incident_all[layer][l];
    }

}



template <typename V, typename L>
void
DynamicInterlayerEdgeStore<V,L>::
erase(
    const L* layer
)
{
    core::assert_not_null(layer, "erase", "layer");

    std::vector<const L*> layers;

    for (auto&& p: edges_)
    {
        layers.push_back(p.first);
    }

    for (auto l: layers)
    {
        edges_[l].erase(layer);
        edge_directionality[l].erase(layer);
        sidx_neighbors_out[l].erase(layer);
        sidx_neighbors_in[l].erase(layer);
        sidx_neighbors_all[l].erase(layer);
        sidx_incident_out[l].erase(layer);
        sidx_incident_in[l].erase(layer);
        sidx_incident_all[l].erase(layer);
    }

    edges_.erase(layer);
    edge_directionality.erase(layer);
    sidx_neighbors_out.erase(layer);
    sidx_neighbors_in.erase(layer);
    sidx_neighbors_all.erase(layer);
    sidx_incident_out.erase(layer);
    sidx_incident_in.erase(layer);
    sidx_incident_all.erase(layer);

}



template <typename V, typename L>
void
DynamicInterlayerEdgeStore<V,L>::
erase(
    const L* layer,
    const V* vertex
)
{
    core::assert_not_null(layer, "erase", "layer");
    core::assert_not_null(vertex, "erase", "vertex");

    std::vector<const L*> layers;

    for (auto&& p: edges_)
    {
        layers.push_back(p.first);
    }

    for (auto l: layers)
    {
        edges_[layer][l];
        edge_directionality[layer].erase(l);
        sidx_neighbors_out[layer].erase(l);
        sidx_neighbors_in[layer].erase(l);
        sidx_neighbors_all[layer].erase(l);
        sidx_incident_out[layer].erase(l);
        sidx_incident_in[layer].erase(l);
        sidx_incident_all[layer].erase(l);
    }

    edges_.erase(layer);
    edge_directionality.erase(layer);
    sidx_neighbors_out.erase(layer);
    sidx_neighbors_in.erase(layer);
    sidx_neighbors_all.erase(layer);
    sidx_incident_out.erase(layer);
    sidx_incident_in.erase(layer);
    sidx_incident_all.erase(layer);

}
}
}

#endif
