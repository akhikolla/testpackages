/**
 * This module defines a generic network co-evolution process.
 *
 * The function "evolve" takes a multilayer network as input and at every step
 * updates each of its layers taking one of the following actions:
 * 1) no action (the layer remains unchanged - used to set different evolution speeds)
 * 2) internal evolution (the layer evolves according to some internal dynamics, defined by an EvolutionModel)
 * 3) external evolution (the layer imports nodes and edges from another layer)
 */

#ifndef UU_GENERATION_PAEVOLUTIONMODEL_H_
#define UU_GENERATION_PAEVOLUTIONMODEL_H_

#include "generation/EvolutionModel.hpp"

namespace uu {
namespace net {


/**
 * @brief Grows a network by first creating a complete graph with m0 nodes, then adding a new node at a time and connecting it to m other nodes chosen with a probability proportional to their degree.
 **/
template <typename M>
class
    PAEvolutionModel :
    public EvolutionModel<M>
{

    size_t m0;
    size_t m;

  public:

    PAEvolutionModel(
        size_t m0,
        size_t m
    );

    ~PAEvolutionModel();

    void
    init_step(
        M* mnet,
        typename M::layer_type* layer,
        GenericObjectList<Vertex>& available_actors
    );

    void
    internal_evolution_step(
        M* mnet,
        typename M::layer_type* layer,
        GenericObjectList<Vertex>& available_actors
    );

    void
    external_evolution_step(
        M* mnet,
        typename M::layer_type* target_layer,
        GenericObjectList<Vertex>& available_actors,
        const typename M::layer_type* ext_layer
    );
};



template <typename M>
PAEvolutionModel<M>::
PAEvolutionModel(
    size_t m0,
    size_t m
)
{
    PAEvolutionModel<M>::m0 = m0;
    PAEvolutionModel<M>::m = m;

    if (m>m0)
    {
        throw core::WrongParameterException("m0 cannot be smaller than m");
    }
}

template <typename M>
PAEvolutionModel<M>::
~PAEvolutionModel(
)
{
}

template <typename M>
void
PAEvolutionModel<M>::
init_step(
    M* mnet,
    typename M::layer_type* layer,
    GenericObjectList<Vertex>& available_actors
)
{
    if (available_actors.size()<m0)
    {
        throw core::WrongParameterException("not enough actors available to initialize the layer (less than m0)");
    }

    std::set<const Vertex*> new_actors;

    // choosing the m0 initial actors
    for (size_t i=0; i<m0; i++)
    {
        auto new_actor = available_actors.get_at_random();
        new_actors.insert(new_actor);
        available_actors.erase(new_actor);
    }

    // adding the actors to the layer
    for (auto actor: new_actors)
    {
        layer->vertices()->add(actor);
    }

    // adding all edges among them
    for (auto node1: *layer->vertices())
    {
        for (auto node2: *layer->vertices())
        {
            if (node1!=node2)
            {
                layer->edges()->add(node1,node2);
            }
        }
    }
}


template <typename M>
void
PAEvolutionModel<M>::
internal_evolution_step(
    M* mnet,
    typename M::layer_type* layer,
    GenericObjectList<Vertex>& available_actors
)
{
    // Choose a new actor to join the layer
    if (available_actors.size()<1)
    {
        return;
    }

    auto new_actor = available_actors.get_at_random();
    available_actors.erase(new_actor);

    layer->vertices()->add(new_actor);

    // Randomly select m nodes with probability proportional to their degree...
    // NOTE: this operation may not be very efficient - think of an alternative implementation
    std::set<const Vertex*> nodes;

    while (nodes.size()<m)
    {
        auto edge = layer->edges()->get_at_random();
        nodes.insert(core::test(.5)?edge->v1:edge->v2);
    }

    // and connect them to the new node
    for (auto old_node: nodes)
    {
        layer->edges()->add(new_actor,old_node);
    }


}


template <typename M>
void
PAEvolutionModel<M>::
external_evolution_step(
    M* mnet,
    typename M::layer_type* target_layer,
    GenericObjectList<Vertex>& available_actors,
    const typename M::layer_type* ext_layer
)
{

    // Choose a new actor to join the layer
    if (available_actors.size()<1)
    {
        return;
    }

    auto new_actor = available_actors.get_at_random();
    available_actors.erase(new_actor);

    target_layer->vertices()->add(new_actor);

    // If the actor is not present on the layer from where we should import contacts, then we do not add any neighbors
    if (!ext_layer->vertices()->contains(new_actor))
    {
        return;
    }

    // now we insert the actor's (at most n) neighbors on ext_layer into the target layer, if they are also present there
    for (auto neighbor: *ext_layer->edges()->neighbors(new_actor,EdgeMode::OUT))
    {
        if (target_layer->vertices()->contains(neighbor))
        {
            target_layer->edges()->add(new_actor,neighbor);
        }
    }
}

}
}

#endif
