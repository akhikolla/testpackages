/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_DATASTRUCTURES_STORES_MDEDGESTORE_H_
#define UU_NET_DATASTRUCTURES_STORES_MDEDGESTORE_H_

#include <memory>
#include <unordered_map>
#include <unordered_map>
#include "core/datastructures/containers/SharedPtrSortedRandomSet.hpp"
#include "core/datastructures/observers/Subject.hpp"
#include "objects/Vertex.hpp"
#include "objects/MLEdge.hpp"
#include "objects/EdgeDir.hpp"
#include "objects/EdgeMode.hpp"
#include "networks/_impl/containers/GenericObjectList.hpp"
#include "core/exceptions/OperationNotSupportedException.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace net {

template <typename VStore>
class MDEdgeStore:
    public core::SharedPtrSortedRandomSet<const MLEdge<Vertex,VStore>>,
            public core::Subject<const MLEdge<Vertex,VStore>>
{

  private:

    typedef core::SharedPtrSortedRandomSet<const MLEdge<Vertex,VStore>> super;

  protected:

    const VStore* layer1;
    const VStore* layer2;

  public:

    typedef MLEdge<Vertex,VStore> value_type;

    /**
     * Constructor.
     */

    MDEdgeStore(
        const VStore* layer1,
        const VStore* layer2,
        EdgeDir dir
    );

  public:

    using super::size;
    using super::add;
    using super::erase;

    virtual
    const MLEdge<Vertex,VStore>*
    add(
        std::shared_ptr<const MLEdge<Vertex,VStore>> e
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
    const MLEdge<Vertex,VStore> *
    add(
        const Vertex* vertex1,
        const VStore* layer1,
        const Vertex* vertex2,
        const VStore* layer2
    );

    virtual
    bool
    erase(
        const MLEdge<Vertex,VStore>* e
    ) override = 0;

    /*
        virtual
        GenericObjectList<MLEdge<Vertex,VStore>>*
                                            get(
                                                const VStore* layer1,
                                                const VStore* layer2
                                            ) const;
    */

    /**
     * @brief Returns the nodes with an edge from/to the input vertex.
     * @param node pointer to the node.
     * @param mode IN, OUT or INOUT.
     * @return the list of neighbors.
     **/
    const
    GenericObjectList<Vertex>*
    neighbors(
        const Vertex* vertex,
        const VStore* layer,
        EdgeMode mode
    ) const;

    /**
     * @brief Returns the nodes with an edge from/to the input vertex.
     * @param node pointer to the node.
     * @param mode IN, OUT or INOUT.
     * @return the list of neighbors.
     **/
    const
    GenericObjectList<MLEdge<Vertex,VStore>>*
                                          incident(
                                                  const Vertex* vertex,
                                                  const VStore* layer,
                                                  EdgeMode mode
                                          ) const;


    bool
    is_directed(
    ) const;

    /*
        void
        set_directed(
            const VStore* layer1,
            const VStore* layer2,
            bool directed
        );
    */

    virtual
    void
    erase(
        const Vertex* vertex,
        const VStore* layer
    );

  protected:


    /** Edges */
    std::unique_ptr<GenericObjectList<MLEdge<Vertex,VStore>>> edges_;

    /** Edge directionality */
    EdgeDir edge_directionality;

    // Indexes to sets of objects (Set IDX):
    std::unordered_map<const VStore*, std::unordered_map<const VStore*, std::unordered_map<const Vertex*, std::unique_ptr<GenericObjectList<Vertex>>>>> sidx_neighbors_out;
    std::unordered_map<const VStore*, std::unordered_map<const VStore*, std::unordered_map<const Vertex*, std::unique_ptr<GenericObjectList<Vertex>>>>> sidx_neighbors_in;
    std::unordered_map<const VStore*, std::unordered_map<const VStore*, std::unordered_map<const Vertex*, std::unique_ptr<GenericObjectList<Vertex>>>>> sidx_neighbors_all;

    std::unordered_map<const VStore*, std::unordered_map<const VStore*, std::unordered_map<const Vertex*, std::unique_ptr<GenericObjectList<MLEdge<Vertex,VStore>>>>>> sidx_incident_out;
    std::unordered_map<const VStore*, std::unordered_map<const VStore*, std::unordered_map<const Vertex*, std::unique_ptr<GenericObjectList<MLEdge<Vertex,VStore>>>>>> sidx_incident_in;
    std::unordered_map<const VStore*, std::unordered_map<const VStore*, std::unordered_map<const Vertex*, std::unique_ptr<GenericObjectList<MLEdge<Vertex,VStore>>>>>> sidx_incident_all;
};


template <typename VStore>
MDEdgeStore<VStore>::
MDEdgeStore(
    const VStore* layer1,
    const VStore* layer2,
    EdgeDir dir
) : layer1(layer1), layer2(layer2), edge_directionality(dir)
{

    core::assert_not_null(layer1, "MDEdgeStore", "layer1");
    core::assert_not_null(layer2, "MDEdgeStore", "layer2");

    edges_ = std::make_unique<GenericObjectList<MLEdge<Vertex,VStore>>>();

    edge_directionality = dir;

    sidx_neighbors_out[layer1][layer2];
    sidx_neighbors_in[layer1][layer2];
    sidx_neighbors_all[layer1][layer2];
    sidx_incident_out[layer1][layer2];
    sidx_incident_in[layer1][layer2];
    sidx_incident_all[layer1][layer2];

    sidx_neighbors_out[layer2][layer1];
    sidx_neighbors_in[layer2][layer1];
    sidx_neighbors_all[layer2][layer1];
    sidx_incident_out[layer2][layer1];
    sidx_incident_in[layer2][layer1];
    sidx_incident_all[layer2][layer1];

}


template <typename VStore>
const MLEdge<Vertex,VStore> *
MDEdgeStore<VStore>::
add(
    const Vertex* vertex1,
    const VStore* layer1,
    const Vertex* vertex2,
    const VStore* layer2
)
{
    core::assert_not_null(vertex1, "add", "vertex1");
    core::assert_not_null(layer1, "add", "layer1");
    core::assert_not_null(vertex2, "add", "vertex2");
    core::assert_not_null(layer2, "add", "layer2");

    auto edge = MLEdge<Vertex,VStore>::create(vertex1, layer1, vertex2, layer2, edge_directionality);
    return add(edge);
}



template <typename VStore>
const MLEdge<Vertex,VStore>*
MDEdgeStore<VStore>::
add(
    std::shared_ptr<const MLEdge<Vertex,VStore>> e
)
{
    core::assert_not_null(e.get(), "add", "e");

    if (e->dir != edge_directionality)
    {
        throw core::OperationNotSupportedException("wrong edge directionality");
    }

    const MLEdge<Vertex,VStore>* new_edge = super::add(e);

    if (!new_edge) // edge already existing
    {
        return nullptr;
    }

    edges_->add(new_edge);

    if (sidx_neighbors_out[e->l1][e->l2].count(e->v1)==0)
    {
        sidx_neighbors_out[e->l1][e->l2][e->v1] = std::make_unique<GenericObjectList<Vertex>>();
        sidx_incident_out[e->l1][e->l2][e->v1] = std::make_unique<GenericObjectList<MLEdge<Vertex,VStore>>>();
    }

    sidx_neighbors_out[e->l1][e->l2][e->v1]->add(e->v2);
    sidx_incident_out[e->l1][e->l2][e->v1]->add(new_edge);


    if (sidx_neighbors_in[e->l2][e->l1].count(e->v2)==0)
    {
        sidx_neighbors_in[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<Vertex>>();
        sidx_incident_in[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<MLEdge<Vertex,VStore>>>();
    }

    sidx_neighbors_in[e->l2][e->l1][e->v2]->add(e->v1);
    sidx_incident_in[e->l2][e->l1][e->v2]->add(new_edge);


    if (sidx_neighbors_all[e->l1][e->l2].count(e->v1)==0)
    {
        sidx_neighbors_all[e->l1][e->l2][e->v1] = std::make_unique<GenericObjectList<Vertex>>();
        sidx_incident_all[e->l1][e->l2][e->v1] = std::make_unique<GenericObjectList<MLEdge<Vertex,VStore>>>();
    }

    sidx_neighbors_all[e->l1][e->l2][e->v1]->add(e->v2);
    sidx_incident_all[e->l1][e->l2][e->v1]->add(new_edge);

    if (sidx_neighbors_all[e->l2][e->l1].count(e->v2)==0)
    {
        sidx_neighbors_all[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<Vertex>>();
        sidx_incident_all[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<MLEdge<Vertex,VStore>>>();
    }

    sidx_neighbors_all[e->l2][e->l1][e->v2]->add(e->v1);
    sidx_incident_all[e->l2][e->l1][e->v2]->add(new_edge);


    if (e->dir == EdgeDir::UNDIRECTED)
    {

        if (sidx_neighbors_out[e->l2][e->l1].count(e->v2)==0)
        {
            sidx_neighbors_out[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<Vertex>>();
            sidx_incident_out[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<MLEdge<Vertex,VStore>>>();
        }

        sidx_neighbors_out[e->l2][e->l1][e->v2]->add(e->v1);
        sidx_incident_out[e->l2][e->l1][e->v2]->add(new_edge);

        if (sidx_neighbors_in[e->l1][e->l2].count(e->v1)==0)
        {
            sidx_neighbors_in[e->l1][e->l2][e->v1] = std::make_unique<GenericObjectList<Vertex>>();
            sidx_incident_in[e->l1][e->l2][e->v1] = std::make_unique<GenericObjectList<MLEdge<Vertex,VStore>>>();
        }

        sidx_neighbors_in[e->l1][e->l2][e->v1]->add(e->v2);
        sidx_incident_in[e->l1][e->l2][e->v1]->add(new_edge);

        /*
        if (sidx_neighbors_all[e->l2][e->l1].count(e->v2)==0)
        {
            sidx_neighbors_all[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<Vertex>>();
            sidx_incident_all[e->l2][e->l1][e->v2] = std::make_unique<GenericObjectList<MLEdge<Vertex,VStore>>>();
        }

        sidx_neighbors_all[e->l2][e->l1][e->v2]->add(e->v1);
        sidx_incident_all[e->l2][e->l1][e->v2]->add(new_edge);
         */
    }

    return new_edge;
}

/*
GenericObjectList<MLEdge<Vertex,VStore>>*
                                    MDEdgeStore<VStore>::
                                    get(
                                        const VStore* layer1,
                                        const VStore* layer2
                                    ) const
{
    core::assert_not_null(layer1, "neighbors", "layer1");
    core::assert_not_null(layer2, "neighbors", "layer2");
    return edges_.get();
}
*/

template <typename VStore>
const
GenericObjectList<Vertex>*
MDEdgeStore<VStore>::
neighbors(
    const Vertex* vertex,
    const VStore* layer,
    EdgeMode mode
) const
{

    core::assert_not_null(layer, "neighbors", "layer");
    core::assert_not_null(vertex, "neighbors", "vertex");

    if (mode==EdgeMode::IN)
    {
        if (sidx_neighbors_in.at(layer).begin()->second.count(vertex)==0)
        {
            return GenericObjectList<Vertex>::empty.get();
        }

        return sidx_neighbors_in.at(layer).begin()->second.at(vertex).get();
    }

    else if (mode==EdgeMode::OUT)
    {
        if (sidx_neighbors_out.at(layer).begin()->second.count(vertex)==0)
        {
            return GenericObjectList<Vertex>::empty.get();
        }

        return sidx_neighbors_out.at(layer).begin()->second.at(vertex).get();
    }

    else if (mode==EdgeMode::INOUT)
    {
        if (sidx_neighbors_all.at(layer).begin()->second.count(vertex)==0)
        {
            return GenericObjectList<Vertex>::empty.get();
        }

        return sidx_neighbors_all.at(layer).begin()->second.at(vertex).get();
    }

    else
    {
        throw core::WrongParameterException("neighborhood mode");
    }
}



template <typename VStore>
const
GenericObjectList<MLEdge<Vertex,VStore>>*
                                      MDEdgeStore<VStore>::
                                      incident(
                                              const Vertex* vertex,
                                              const VStore* layer,
                                              EdgeMode mode
                                      ) const
{

    core::assert_not_null(layer1, "neighbors", "layer");
    core::assert_not_null(vertex, "incident", "vertex");

    if (mode==EdgeMode::IN)
    {
        if (sidx_incident_in.at(layer).begin()->second.count(vertex)==0)
        {
            return GenericObjectList<MLEdge<Vertex,VStore>>::empty.get();
        }

        return sidx_incident_in.at(layer).begin()->second.at(vertex).get();
    }

    else if (mode==EdgeMode::OUT)
    {
        if (sidx_incident_out.at(layer).begin()->second.count(vertex)==0)
        {
            return GenericObjectList<MLEdge<Vertex,VStore>>::empty.get();
        }

        return sidx_incident_out.at(layer).begin()->second.at(vertex).get();
    }

    else if (mode==EdgeMode::INOUT)
    {
        if (sidx_incident_all.at(layer).begin()->second.count(vertex)==0)
        {
            return GenericObjectList<MLEdge<Vertex,VStore>>::empty.get();
        }

        return sidx_incident_all.at(layer).begin()->second.at(vertex).get();
    }

    else
    {
        throw core::WrongParameterException("neighborhood mode");
    }
}



template <typename VStore>
bool
MDEdgeStore<VStore>::
is_directed(
) const
{
    return edge_directionality == EdgeDir::DIRECTED?true:false;
}



/*
void
MDEdgeStore<VStore>::
set_directed(
const VStore* layer1,
const VStore* layer2,
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
*/

// @todo What?...

template <typename VStore>
void
MDEdgeStore<VStore>::
erase(
    const Vertex* vertex,
    const VStore* layer
)
{
    core::assert_not_null(layer, "erase", "layer");
    core::assert_not_null(vertex, "erase", "vertex");

    std::unordered_set<const MLEdge<Vertex,VStore>*> to_erase;

    for (auto e: *incident(vertex,layer,EdgeMode::INOUT))
    {
        to_erase.insert(e);
    }

    for (auto e: to_erase)
    {
        erase(e);
    }

}

}
}

#endif
