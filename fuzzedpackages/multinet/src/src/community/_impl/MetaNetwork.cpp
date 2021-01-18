
#include "community/_impl/MetaNetwork.hpp"
#include <unordered_map>
#include <vector>
#include <memory>
#include "networks/WeightedNetwork.hpp"
#include "objects/EdgeMode.hpp"
#include "objects/Vertex.hpp"
#include "objects/Edge.hpp"

namespace uu {
namespace net {


MetaNetwork::
MetaNetwork()
{
    w = std::make_unique<WeightedNetwork>("w", EdgeDir::UNDIRECTED, true);
}

const Vertex*
MetaNetwork::
add(
    const Vertex* u
)
{
    auto v = w->vertices()->add(std::to_string(order));
    order++;
    mapping[v];
    mapping[v].insert(u);
    reverse_mapping[u] = v;
    return v;
}


const Edge*
MetaNetwork::
edge(
    const Vertex* u,
    const Vertex* v,
    double weight
)
{
    auto u_prime = reverse_mapping.at(u);
    auto v_prime = reverse_mapping.at(v);

    auto e = w->edges()->add(u_prime, v_prime);

    double previous_weight = 0.0;

    if (!e)
    {
        e = w->edges()->get(u_prime, v_prime);
        previous_weight = w->get_weight(e).value;
    }

    w->set_weight(e, previous_weight+weight);
    return e;
}

const WeightedNetwork*
MetaNetwork::
get(
) const
{
    return w.get();
}



}
}
