#ifndef UU_MNET_COMMUNITY_PILLARCOMMUNITY_H_
#define UU_MNET_COMMUNITY_PILLARCOMMUNITY_H_

#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <utility>
#include "objects/EdgeMode.hpp"
#include "community/VertexLayerCommunity.hpp"

namespace uu {
namespace net {

template <typename L>
class
    PillarCommunity
{
  public:

    std::string
    to_string(
    ) const;

    void
    add_actor(
        const Vertex*
    );

    const std::unordered_set<const Vertex*>&
    get_actors(
    ) const;

    size_t
    num_actors(
    ) const;

    void
    add_layer(
        const L*
    );

    const std::unordered_set<const L*>&
    get_layers(
    ) const;

    size_t
    num_layers(
    ) const;

  private:
    std::unordered_set<const Vertex*> actors;
    std::unordered_set<const L*> layers;
};

template <typename L>
std::unique_ptr<VertexLayerCommunity<const L>>
        to_vertex_layer_community(
            PillarCommunity<L>* com
        );


template <typename L>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const L>>>
to_vertex_layer_community_structure(
    CommunityStructure<PillarCommunity<L>>* com
);
// Definitions


template <typename L>
std::unique_ptr<VertexLayerCommunity<const L>>
        to_vertex_layer_community(
            PillarCommunity<L>* com
        )
{
    auto res = std::make_unique<VertexLayerCommunity<const L>>();

    for (auto a: com->get_actors())
    {
        for (auto l: com->get_layers())
        {
            res->add(std::make_pair(a,l));
        }
    }

    return res;
}

template <typename L>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const L>>>
to_vertex_layer_community_structure(
    CommunityStructure<PillarCommunity<L>>* com
)
{
    auto res = std::make_unique<CommunityStructure<VertexLayerCommunity<const L>>>();

    for (auto c: *com)
    {
        res->add(std::move(to_vertex_layer_community(c)));
    }

    return res;
}

template <typename L>
std::string
PillarCommunity<L>::
to_string(
) const
{
    std::string result = "";
    size_t idx = 0;

    for (auto actor: actors)
    {
        if (idx==0)
        {
            result += "[";
        }

        result += actor->to_string();

        if (idx!=actors.size()-1)
        {
            result += ", ";
        }

        else
        {
            result += "]";
        }

        idx++;
    }

    for (auto layer: layers)
    {
        if (idx==0)
        {
            result += "::[";
        }

        result += layer->to_string();

        if (idx!=layers.size()-1)
        {
            result += ", ";
        }

        else
        {
            result += "]";
        }

        idx++;
    }

    return result;
}


template <typename L>
void PillarCommunity<L>::
add_actor(
    const Vertex* actor
)
{
    actors.insert(actor);
}


template <typename L>
const std::unordered_set<const Vertex*>&
PillarCommunity<L>::
get_actors(
) const
{
    return actors;
}


template <typename L>
void
PillarCommunity<L>::
add_layer(
    const L* layer
)
{
    layers.insert(layer);
}


template <typename L>
const std::unordered_set<const L*>&
PillarCommunity<L>::
get_layers(
) const
{
    return layers;
}

}
}


#endif
