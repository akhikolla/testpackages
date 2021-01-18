#ifndef UU_MNET_COMMUNITY_MLCPMCOMMUNITY_H_
#define UU_MNET_COMMUNITY_MLCPMCOMMUNITY_H_

#include "objects/MultiplexClique.hpp"
#include "community/VertexLayerCommunity.hpp"


namespace uu {
namespace net {

// DATA STRUCTURES

template <typename M>
struct layer_set_comparator
{
    bool
    operator()(
        const typename std::set<const typename M::layer_type*>& a,
        const typename std::set<const typename M::layer_type*>& b
    ) const
    {
        if (a.size() != b.size())
        {
            return a.size() < b.size();
        }

        typename std::set<const typename M::layer_type*>::iterator it1 = a.begin();
        typename std::set<const typename M::layer_type*>::iterator it2 = b.begin();

        for (size_t i=0; i<a.size(); i++)
        {
            if ((*it1)<(*it2))
            {
                return true;
            }

            if ((*it1)>(*it2))
            {
                return false;
            }

            ++it1;
            ++it2;
        }

        return false;
    }
};

template <typename M>
using layer_sets = std::set<typename std::set<const typename M::layer_type*>,layer_set_comparator<M> >;


/* the ML-CPM algorithm uses a special type of community, defined as a set of adjacent cliques */
template <typename M>
class
    MLCPMCommunity
{
  public:
    MLCPMCommunity(
    );

    MLCPMCommunity(
        size_t cid,
        std::unordered_set<std::shared_ptr<MultiplexClique<M>>> cliques,
        typename std::unordered_set<const typename M::layer_type*> layers
    );

    static
    std::shared_ptr<MLCPMCommunity<M>>
                                    create(
                                    );

    void
    add_clique(
        std::shared_ptr<MultiplexClique<M>>
    );

    void
    add_layer(
        const typename M::layer_type*
    );

    const
    typename std::set<const typename M::layer_type*>&
    get_layers(
    );

    std::set<const Vertex*>
    actors(
    ) const;

    std::unique_ptr<VertexLayerCommunity<const typename M::layer_type>>
            to_community(
                //const M* net
            ) const;

    int
    size(
    ) const;

    bool
    operator==(
        const MLCPMCommunity& comp
    ) const;

    bool
    operator!=(
        const MLCPMCommunity& comp
    ) const;

    bool
    operator<(
        const MLCPMCommunity& comp
    ) const;

    bool
    operator>(
        const MLCPMCommunity& comp
    ) const;

    std::string
    to_string();

    long id;

    std::set<std::shared_ptr<MultiplexClique<M>>> cliques;

    typename std::set<const typename M::layer_type*> layers;

};


template <typename M>
MLCPMCommunity<M>::
MLCPMCommunity() :
    id(0)
{

}

template <typename M>
MLCPMCommunity<M>::
MLCPMCommunity(
    size_t id,
    std::unordered_set<std::shared_ptr<MultiplexClique<M>>> cliques,
    std::unordered_set<const typename M::layer_type*> layers
) :
    id(id), cliques(cliques.begin(),cliques.end()), layers(layers.begin(),layers.end()) {}

template <typename M>
std::set<const Vertex*>
MLCPMCommunity<M>::
actors(
) const
{
    std::set<const Vertex*> actors;

    for (std::shared_ptr<MultiplexClique<M>> clique: cliques)
    {
        for (auto actor: clique->actors)
        {
            actors.insert(actor);
        }
    }

    return actors;
}

template <typename M>
std::shared_ptr<MLCPMCommunity<M>>
                                MLCPMCommunity<M>::
                                create(
                                )
{
    return std::make_shared<MLCPMCommunity<M>>();
}

template <typename M>
void
MLCPMCommunity<M>::
add_clique(std::shared_ptr<MultiplexClique<M>> clique)
{
    cliques.insert(clique);
}

template <typename M>
void
MLCPMCommunity<M>::
add_layer(
    const typename M::layer_type* layer
)
{
    layers.insert(layer);
}

template <typename M>
const
typename std::set<const typename M::layer_type*>&
MLCPMCommunity<M>::
get_layers(
)
{
    return layers;
}

template <typename M>
std::unique_ptr<VertexLayerCommunity<const typename M::layer_type>>
        MLCPMCommunity<M>::
        to_community(
            //const M* net
        ) const
{
    auto result = std::make_unique<VertexLayerCommunity<const typename M::layer_type>>();

    for (auto actor: actors())
    {
        for (auto layer: layers)
        {
            if (!layer->vertices()->contains(actor))
            {
                continue;
            }

            auto iv = std::make_pair(actor, layer);
            result->add(iv);
        }
    }

    return result;
}

template <typename M>
std::string
MLCPMCommunity<M>::
to_string(
)
{
    std::string res = "C" + std::to_string(id) + ": ";
    std::unordered_set<const Vertex*> actors;

    for (std::shared_ptr<MultiplexClique<M>> clique: cliques)
    {
        for (auto actor: clique->actors)
        {
            actors.insert(actor);
        }
    }

    for (auto actor: actors)
    {
        res += actor->name + " ";
    }

    res += "( ";

    for (const typename M::layer_type* layer: layers)
    {
        res += layer->name + " ";
    }

    res += ")";
    return res;
}

template <typename M>
int
MLCPMCommunity<M>::
size(
) const
{
    std::unordered_set<const Vertex*> actors;

    for (std::shared_ptr<MultiplexClique<M>> clique: cliques)
    {
        for (auto actor: clique->actors)
        {
            actors.insert(actor);
        }
    }

    return actors.size();
}

}
}


#endif
