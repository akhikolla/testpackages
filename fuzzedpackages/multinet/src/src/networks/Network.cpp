#include "networks/Network.hpp"

#include "networks/_impl/observers/NoLoopCheckObserver.hpp"

namespace uu {
namespace net {

Network::
Network(
    const std::string& name,
    EdgeDir dir,
    bool allows_loops) : name(name)
{

    auto vs = std::make_unique<AttrVertexStore>();

    auto es = std::make_unique<AttrSimpleEdgeStore>(dir);

    GraphType t;
    t.allows_loops = allows_loops;
    t.is_directed = dir==EdgeDir::DIRECTED ? true : false;
    t.is_weighted = false;

    data_ = std::make_unique<Graph<AttrVertexStore, AttrSimpleEdgeStore>>(name, t, std::move(vs), std::move(es));

    if (!allows_loops)
    {
        auto obs = std::make_unique<NoLoopCheckObserver>();
        data_->edges()->attach(obs.get());
        data_->register_observer(std::move(obs));
    }
}



AttrVertexStore*
Network::
vertices(
)
{
    return data_->vertices();
}



const AttrVertexStore*
Network::
vertices(
) const
{
    return data_->vertices();
}


AttrSimpleEdgeStore*
Network::
edges(
)
{
    return data_->edges();
}


const AttrSimpleEdgeStore*
Network::
edges(
) const
{
    return data_->edges();
}


bool
Network::
is_directed(
) const
{
    return data_->is_directed();
}


bool
Network::
is_weighted(
) const
{
    return data_->is_weighted();
}


bool
Network::
is_probabilistic(
) const
{
    return data_->is_probabilistic();
}


bool
Network::
is_temporal(
) const
{
    return data_->is_temporal();
}


bool
Network::
is_attributed(
) const
{
    return data_->is_attributed();
}


bool
Network::
allows_multi_edges(
) const
{
    return data_->allows_multi_edges();
}


bool
Network::
allows_loops(
) const
{
    return data_->allows_loops();
}

}
}

