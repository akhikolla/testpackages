/**

 */

#include "networks/_impl/stores/MultiEdgeStore.hpp"

namespace uu {
namespace net {


MultiEdgeStore::
MultiEdgeStore(
    EdgeDir dir
) :
    EdgeStore(
        dir
    )
{
}


const Edge *
MultiEdgeStore::
add(
    const Vertex* vertex1,
    const Vertex* vertex2
)
{

    std::shared_ptr<const Edge> e = Edge::create(vertex1, vertex2, edge_directionality);

    return add(e);
}


const Edge*
MultiEdgeStore::
add(
    std::shared_ptr<const Edge> e
)
{
    for (auto obs: observers)
    {
        obs->notify_add(e.get());
    }

    // No need to check for edge existence

    // EDGE CREATION

    const Edge * new_edge = super::add(e);

    if (!new_edge)
    {
        return nullptr;
    }


    /// MULTI SPEC.
    cidx_edges_by_vertices[e->v1][e->v2].insert(new_edge);

    /// DIR SPECIFIC.

    if (!is_directed())
    {
        cidx_edges_by_vertices[e->v2][e->v1].insert(new_edge);

    }

    return new_edge;
}



core::SortedRandomSet<const Edge*>
MultiEdgeStore::
get(
    const Vertex* vertex1,
    const Vertex* vertex2
) const
{
    core::SortedRandomSet<const Edge*> result;

    if (cidx_edges_by_vertices.count(vertex1)>0 &&
            cidx_edges_by_vertices.at(vertex1).count(vertex2)>0)
    {
        auto edges = cidx_edges_by_vertices.at(vertex1).at(vertex2);

        for (auto edge: edges)
        {
            result.add(edge);
        }
    }

    return result;

}





bool
MultiEdgeStore::
erase(
    const Edge* edge
)
{
    for (auto obs: observers)
    {
        obs->notify_erase(edge);
    }

    cidx_edges_by_vertices[edge->v1][edge->v2].erase(edge);

    if (cidx_edges_by_vertices[edge->v1][edge->v2].size()==0)
    {
        sidx_neighbors_in[edge->v2]->erase(edge->v1);
        sidx_neighbors_out[edge->v1]->erase(edge->v2);
        sidx_incident_in[edge->v2]->erase(edge);
        sidx_incident_out[edge->v1]->erase(edge);
    }

    // if the edge is directed, we erase neighbors only if there isn't
    // any edge in the other direction keeping them neighbors
    if (edge->dir==EdgeDir::DIRECTED && cidx_edges_by_vertices[edge->v2][edge->v1].size()==0)
    {
        sidx_neighbors_all[edge->v2]->erase(edge->v1);
        sidx_neighbors_all[edge->v1]->erase(edge->v2);
        sidx_incident_all[edge->v2]->erase(edge);
        sidx_incident_all[edge->v1]->erase(edge);
    }

    if (edge->dir==EdgeDir::UNDIRECTED)
    {
        cidx_edges_by_vertices[edge->v2][edge->v1].erase(edge);

        if (cidx_edges_by_vertices[edge->v1][edge->v2].size()==0)
        {
            sidx_neighbors_in[edge->v1]->erase(edge->v2);
            sidx_neighbors_out[edge->v2]->erase(edge->v1);
            sidx_neighbors_all[edge->v1]->erase(edge->v2);
            sidx_neighbors_all[edge->v2]->erase(edge->v1);
            sidx_incident_in[edge->v1]->erase(edge);
            sidx_incident_out[edge->v2]->erase(edge);
            sidx_incident_all[edge->v1]->erase(edge);
            sidx_incident_all[edge->v2]->erase(edge);
        }
    }

    bool res = core::SharedPtrSortedRandomSet<const Edge>::erase(edge);
    return res;
}


void
MultiEdgeStore::
erase(
    const Vertex* vertex
)
{
    std::unordered_set<const Edge*> to_erase;

    for (const Vertex* neighbor: *neighbors(vertex,EdgeMode::OUT))
    {
        auto edges = get(vertex,neighbor);

        for (auto edge: edges)
        {
            to_erase.insert(edge);
        }
    }

    for (const Vertex* neighbor: *neighbors(vertex,EdgeMode::IN))
    {
        auto edges = get(neighbor,vertex);

        for (auto edge: edges)
        {
            to_erase.insert(edge);
        }
    }

    for (auto e: to_erase)
    {
        erase(e);
    }
}


std::string
MultiEdgeStore::
summary(
) const
{
    std::string summary = std::to_string(size()) + " multiedges";
    return summary;
}



}
}
