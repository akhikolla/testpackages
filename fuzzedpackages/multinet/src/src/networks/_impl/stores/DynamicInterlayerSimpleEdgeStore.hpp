/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_DATASTRUCTURES_STORES_DYNAMICINTERLAYERSIMPLEEDGESTORE_H_
#define UU_NET_DATASTRUCTURES_STORES_DYNAMICINTERLAYERSIMPLEEDGESTORE_H_

#include "networks/_impl/stores/DynamicInterlayerEdgeStore.hpp"
#include "networks/Network.hpp"
#include "networks/_impl/stores/SimpleEdgeStore.hpp"
#include "core/exceptions/ElementNotFoundException.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace net {

template <typename V, typename L>
class
    DynamicInterlayerSimpleEdgeStore :
    public DynamicInterlayerEdgeStore<V,L>
{
    typedef DynamicInterlayerEdgeStore<V,L> super;

  public:
    using super::super;

    using super::add;
    using super::get;
    using super::neighbors;
    using super::is_directed;
    using super::attach;
    using super::erase;
    using super::size;
    //using super::summary;

    using super::edges_;
    using super::edge_directionality;
    using super::sidx_neighbors_out;
    using super::sidx_neighbors_in;
    using super::sidx_neighbors_all;
    using super::sidx_incident_out;
    using super::sidx_incident_in;
    using super::sidx_incident_all;
    using super::observers;



    /**
     * Adds a new edge.
     * @param e edge to be added.
     * @return a pointer to the new edge, or nullptr if the edge already exists.
     **/
    virtual
    const MLEdge<V,L>*
    add(
        std::shared_ptr<const MLEdge<V,L>>  e
    ) override;


    /**
     * Returns an edge.
     * This function can also be used to check if an edge is present.
     * @param vertex1 a pointer to the "from" actor if directed, or to one end
     * of the edge if undirected.
     * @param vertex2 a pointer to the "to" actor if directed, or one end
     * of the edge if undirected.
     * @return a pointer to the requested edge, or nullptr if it does not exist.
     **/
    const MLEdge<V,L>*
    get(
        const V* vertex1,
        const L* layer1,
        const V* vertex2,
        const L* layer2
    ) const;



    virtual
    bool
    erase(
        const MLEdge<V,L>* e
    ) override;


    virtual
    void
    add(
        const L* layer
    ) override;

    virtual
    void
    erase(
        const L* layer
    ) override;


    virtual
    void
    erase(
        const L* layer,
        const V* vertex
    ) override;


  protected:

    // Indexes to objects (Component IDX):
    std::unordered_map<const L*, std::unordered_map<const L*, std::unordered_map<const V*, std::unordered_map<const V*, const MLEdge<V,L>*>>>> cidx_edge_by_vertexes;
};


template<typename V, typename L>
const MLEdge<V,L> *
DynamicInterlayerSimpleEdgeStore<V,L>::
add(
    std::shared_ptr<const MLEdge<V,L>> e
)
{
    core::assert_not_null(e.get(), "add", "e");

    for (auto obs: observers)
    {
        obs->notify_add(e.get());
    }

    // get() also checks if the layers are present in this store
    if (get(e->v1, e->l1, e->v2, e->l2))
    {
        return nullptr;
    }

    auto new_edge = super::add(e);

    if (!new_edge)
    {
        return nullptr;
    }

    cidx_edge_by_vertexes[e->l1][e->l2][e->v1][e->v2] = new_edge;

    /// DIR SPECIFIC.

    if (!is_directed(e->l1, e->l2))
    {
        cidx_edge_by_vertexes[e->l2][e->l1][e->v2][e->v1] = new_edge;
    }

    return new_edge;
}


template<typename V, typename L>
const MLEdge<V,L>*
DynamicInterlayerSimpleEdgeStore<V,L>::
get(
    const V* vertex1,
    const L* layer1,
    const V* vertex2,
    const L* layer2
) const
{

    core::assert_not_null(vertex1, "get", "vertex1");
    core::assert_not_null(layer1, "get", "layer1");
    core::assert_not_null(vertex2, "get", "vertex2");
    core::assert_not_null(layer2, "get", "layer2");

    auto l1 = cidx_edge_by_vertexes.find(layer1);

    if (l1 == cidx_edge_by_vertexes.end())
    {
        throw core::ElementNotFoundException("layer " + layer1->name +
                                             " is not present in this store");
    }

    auto l2 = l1->second.find(layer2);

    if (l2 == l1->second.end())
    {
        throw core::ElementNotFoundException("layer " + layer2->name +
                                             " is not present in this store");
    }

    auto v1 = l2->second.find(vertex1);

    if (v1 == l2->second.end())
    {
        return nullptr;
    }

    auto v2 = v1->second.find(vertex2);

    if (v2 == v1->second.end())
    {
        return nullptr;
    }

    else
    {
        return v2->second;
    }
}



template<typename V, typename L>
bool
DynamicInterlayerSimpleEdgeStore<V,L>::
erase(
    const MLEdge<V,L>* edge
)
{
    core::assert_not_null(edge, "erase", "edge");

    for (auto obs: observers)
    {
        obs->notify_erase(edge);
    }

    edges_[edge->l1][edge->l2]->erase(edge);
    edges_[edge->l2][edge->l1]->erase(edge);

    cidx_edge_by_vertexes[edge->l1][edge->l2][edge->v1].erase(edge->v2);

    sidx_neighbors_in[edge->l2][edge->l1][edge->v2]->erase(edge->v1);
    sidx_neighbors_out[edge->l1][edge->l2][edge->v1]->erase(edge->v2);
    sidx_incident_in[edge->l2][edge->l1][edge->v2]->erase(edge);
    sidx_incident_out[edge->l1][edge->l2][edge->v1]->erase(edge);


    // if the edge is directed, we erase neighbors only if there isn't
    // an edge in the other direction keeping them neighbors
    if (is_directed(edge->l1, edge->l2))
    {

        if (!get(edge->v2,edge->l2,edge->v1,edge->l1))
        {
            sidx_neighbors_all[edge->l2][edge->l1][edge->v2]->erase(edge->v1);
            sidx_neighbors_all[edge->l1][edge->l2][edge->v1]->erase(edge->v2);
            sidx_incident_all[edge->l2][edge->l1][edge->v2]->erase(edge);
            sidx_incident_all[edge->l1][edge->l2][edge->v1]->erase(edge);
        }
    }

    else
    {

        cidx_edge_by_vertexes[edge->l2][edge->l1][edge->v2].erase(edge->v1);

        sidx_neighbors_in[edge->l1][edge->l2][edge->v1]->erase(edge->v2);
        sidx_neighbors_out[edge->l2][edge->l1][edge->v2]->erase(edge->v1);
        sidx_neighbors_all[edge->l1][edge->l2][edge->v1]->erase(edge->v2);
        sidx_neighbors_all[edge->l2][edge->l1][edge->v2]->erase(edge->v1);
        sidx_incident_in[edge->l1][edge->l2][edge->v1]->erase(edge);
        sidx_incident_out[edge->l2][edge->l1][edge->v2]->erase(edge);
        sidx_incident_all[edge->l1][edge->l2][edge->v1]->erase(edge);
        sidx_incident_all[edge->l2][edge->l1][edge->v2]->erase(edge);
    }


    return core::SharedPtrSortedRandomSet<const MLEdge<V,L>>::erase(edge);
}



template<typename V, typename L>
void
DynamicInterlayerSimpleEdgeStore<V,L>::
erase(
    const L* layer,
    const V* vertex
)
{

    core::assert_not_null(layer, "erase", "layer");
    core::assert_not_null(vertex, "erase", "vertex");

    std::unordered_set<const MLEdge<V,L>*> to_erase;

    std::vector<const L*> layers;

    for (auto&& p: edges_)
    {
        layers.push_back(p.first);
    }



    for (auto l: layers)
    {
        for (auto neighbor: *neighbors(layer,l,vertex,EdgeMode::OUT))
        {

            auto e = get(vertex,layer,neighbor,l);

            to_erase.insert(e);
        }
    }

    for (auto l: layers)
    {
        for (auto neighbor: *neighbors(layer,l,vertex,EdgeMode::IN))
        {
            auto e = get(neighbor,l,vertex,layer);

            to_erase.insert(e);
        }
    }


    for (auto e: to_erase)
    {
        erase(e);
    }
}


template <typename V, typename L>
void
DynamicInterlayerSimpleEdgeStore<V,L>::
add(
    const L* layer
)
{
    super::add(layer);

    std::vector<const L*> layers;

    for (auto&& p: edges_)
    {
        layers.push_back(p.first);
    }

    for (auto l: layers)
    {
        cidx_edge_by_vertexes[l][layer];
        cidx_edge_by_vertexes[layer][l];
    }

}

template <typename V, typename L>
void
DynamicInterlayerSimpleEdgeStore<V,L>::
erase(
    const L* layer
)
{
    super::erase(layer);


    std::vector<const L*> layers;

    for (auto&& p: edges_)
    {
        layers.push_back(p.first);
    }

    for (auto l: layers)
    {
        cidx_edge_by_vertexes[l].erase(layer);
    }

    cidx_edge_by_vertexes.erase(layer);
}


}
}

#endif
