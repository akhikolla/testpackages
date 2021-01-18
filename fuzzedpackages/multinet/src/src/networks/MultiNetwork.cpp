#include "networks/MultiNetwork.hpp"

#include "networks/_impl/observers/NoLoopCheckObserver.hpp"

namespace uu {
namespace net {

MultiNetwork::
MultiNetwork(
    const std::string& name,
    EdgeDir dir,
    bool allows_loops) : name(name)
{

    auto vs = std::make_unique<AttrVertexStore>();

    auto es = std::make_unique<AttrMultiEdgeStore>(dir);

    GraphType t;
    t.allows_loops = allows_loops;
    t.allows_multi_edges = true;
    t.is_directed = dir==EdgeDir::DIRECTED ? true : false;
    t.is_weighted = false;

    data_ = std::make_unique<Graph<AttrVertexStore, AttrMultiEdgeStore>>(name, t, std::move(vs), std::move(es));

    if (!allows_loops)
    {
        auto obs = std::make_unique<NoLoopCheckObserver>();
        data_->edges()->attach(obs.get());
        data_->register_observer(std::move(obs));
    }
}



AttrVertexStore*
MultiNetwork::
vertices(
)
{
    return data_->vertices();
}


const AttrVertexStore*
MultiNetwork::
vertices(
) const
{
    return data_->vertices();
}


AttrMultiEdgeStore*
MultiNetwork::
edges(
)
{
    return data_->edges();
}


const AttrMultiEdgeStore*
MultiNetwork::
edges(
) const
{
    return data_->edges();
}


bool
MultiNetwork::
is_directed(
) const
{
    return data_->is_directed();
}


bool
MultiNetwork::
allows_loops(
) const
{
    return data_->allows_loops();
}


bool
MultiNetwork::
is_weighted(
) const
{
    return data_->is_weighted();
}


bool
MultiNetwork::
is_probabilistic(
) const
{
    return data_->is_probabilistic();
}


bool
MultiNetwork::
is_temporal(
) const
{
    return data_->is_temporal();
}


bool
MultiNetwork::
is_attributed(
) const
{
    return data_->is_attributed();
}


bool
MultiNetwork::
allows_multi_edges(
) const
{
    return data_->allows_multi_edges();
}

}
}

