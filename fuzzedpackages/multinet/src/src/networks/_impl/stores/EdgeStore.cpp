/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#include <unordered_set>
#include <unordered_map>
#include "networks/_impl/stores/EdgeStore.hpp"

namespace uu {
namespace net {


EdgeStore::
EdgeStore(
    EdgeDir dir
)
{
    edge_directionality = dir;
}


const Edge*
EdgeStore::
add(
    std::shared_ptr<const Edge> e
)
{
    core::assert_not_null(e.get(), "add", "e");

    if (e->dir != edge_directionality)
    {
        throw core::OperationNotSupportedException("wrong edge directionality");
    }

    const Edge* new_edge = core::SharedPtrSortedRandomSet<const Edge>::add(e);

    if (!new_edge) // edge already existing
    {
        return nullptr;
    }

    if (sidx_neighbors_out.count(e->v1)==0)
    {
        sidx_neighbors_out[e->v1] = std::make_unique<GenericObjectList<Vertex>>();
        sidx_incident_out[e->v1] = std::make_unique<GenericObjectList<Edge>>();
    }

    sidx_neighbors_out[e->v1]->add(e->v2);
    sidx_incident_out[e->v1]->add(new_edge);

    if (sidx_neighbors_in.count(e->v2)==0)
    {
        sidx_neighbors_in[e->v2] = std::make_unique<GenericObjectList<Vertex>>();
        sidx_incident_in[e->v2] = std::make_unique<GenericObjectList<Edge>>();
    }

    sidx_neighbors_in[e->v2]->add(e->v1);
    sidx_incident_in[e->v2]->add(new_edge);

    if (sidx_neighbors_all.count(e->v1)==0)
    {
        sidx_neighbors_all[e->v1] = std::make_unique<GenericObjectList<Vertex>>();
        sidx_incident_all[e->v1] = std::make_unique<GenericObjectList<Edge>>();
    }

    sidx_neighbors_all[e->v1]->add(e->v2);
    sidx_incident_all[e->v1]->add(new_edge);

    if (sidx_neighbors_all.count(e->v2)==0)
    {
        sidx_neighbors_all[e->v2] = std::make_unique<GenericObjectList<Vertex>>();
        sidx_incident_all[e->v2] = std::make_unique<GenericObjectList<Edge>>();
    }

    sidx_neighbors_all[e->v2]->add(e->v1);
    sidx_incident_all[e->v2]->add(new_edge);

    if (!is_directed())
    {

        if (sidx_neighbors_out.count(e->v2)==0)
        {
            sidx_neighbors_out[e->v2] = std::make_unique<GenericObjectList<Vertex>>();
            sidx_incident_out[e->v2] = std::make_unique<GenericObjectList<Edge>>();
        }

        sidx_neighbors_out[e->v2]->add(e->v1);
        sidx_incident_out[e->v2]->add(new_edge);

        if (sidx_neighbors_in.count(e->v1)==0)
        {
            sidx_neighbors_in[e->v1] = std::make_unique<GenericObjectList<Vertex>>();
            sidx_incident_in[e->v1] = std::make_unique<GenericObjectList<Edge>>();
        }

        sidx_neighbors_in[e->v1]->add(e->v2);
        sidx_incident_in[e->v1]->add(new_edge);
    }

    return new_edge;
}



const
GenericObjectList<Vertex>*
EdgeStore::
neighbors(
    const Vertex* vertex,
    EdgeMode mode
) const
{

    core::assert_not_null(vertex, "neighbors", "vertex");

    if (mode==EdgeMode::IN)
    {
        if (sidx_neighbors_in.count(vertex)==0)
        {
            return GenericObjectList<Vertex>::empty.get();
        }

        return sidx_neighbors_in.at(vertex).get();
    }

    else if (mode==EdgeMode::OUT)
    {
        if (sidx_neighbors_out.count(vertex)==0)
        {
            return GenericObjectList<Vertex>::empty.get();
        }

        return sidx_neighbors_out.at(vertex).get();
    }

    else if (mode==EdgeMode::INOUT)
    {
        if (sidx_neighbors_all.count(vertex)==0)
        {
            return GenericObjectList<Vertex>::empty.get();
        }

        return sidx_neighbors_all.at(vertex).get();
    }

    else
    {
        throw core::WrongParameterException("neighborhood mode");
    }
}


const
GenericObjectList<Edge>*
EdgeStore::
incident(
    const Vertex* vertex,
    EdgeMode mode
) const
{

    core::assert_not_null(vertex, "incident", "vertex");

    if (mode==EdgeMode::IN)
    {
        if (sidx_incident_in.count(vertex)==0)
        {
            return GenericObjectList<Edge>::empty.get();
        }

        return sidx_incident_in.at(vertex).get();
    }

    else if (mode==EdgeMode::OUT)
    {
        if (sidx_incident_out.count(vertex)==0)
        {
            return GenericObjectList<Edge>::empty.get();
        }

        return sidx_incident_out.at(vertex).get();
    }

    else if (mode==EdgeMode::INOUT)
    {
        if (sidx_incident_all.count(vertex)==0)
        {
            return GenericObjectList<Edge>::empty.get();
        }

        return sidx_incident_all.at(vertex).get();
    }

    else
    {
        throw core::WrongParameterException("neighborhood mode");
    }
}


bool
EdgeStore::
is_directed(
)
{
    return edge_directionality==EdgeDir::DIRECTED?true:false;
}


std::string
EdgeStore::
summary(
) const
{
    size_t s = size();

    std::string summary = std::to_string(s) +
                          (s==1?" edge":" edges");
    return summary;
}

}
}
