
#include "community/_impl/GMetaNetwork.hpp"
#include <unordered_map>
#include <vector>
#include <memory>
#include "networks/WeightedNetwork.hpp"
#include "objects/EdgeMode.hpp"
#include "objects/Vertex.hpp"
#include "objects/Edge.hpp"

namespace uu {
namespace net {


GMetaNetwork::
GMetaNetwork()
{
    w = std::make_unique<MultiNetwork>("w", EdgeDir::UNDIRECTED, true);
}

const Vertex*
GMetaNetwork::
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
GMetaNetwork::
edge(
    const Vertex* u,
    const Vertex* v,
    size_t type,
    double weight
)
{
    auto u_prime = reverse_mapping.at(u);
    auto v_prime = reverse_mapping.at(v);

    for (auto edge: w->edges()->get(u_prime, v_prime))
    {
        if (edge_type.at(edge) == type)
        {
            edge_weight[edge] += weight;
            return edge;
        }
    }

    // if no edge of type type exists:

    auto e = w->edges()->add(u_prime, v_prime);
    edge_type[e] = type;
    edge_weight[e] = weight;

    return e;
}

double
GMetaNetwork::
get_weight(
    const Edge* e
) const
{
    return edge_weight.at(e);
}

size_t
GMetaNetwork::
get_type(
    const Edge* e
) const
{
    return edge_type.at(e);
}

const MultiNetwork*
GMetaNetwork::
get(
) const
{
    return w.get();
}



}
}
