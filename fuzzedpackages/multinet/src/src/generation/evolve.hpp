/**
 * This module defines a generic network co-evolution process.
 *
 * The function "evolve" takes a multilayer network as input and at every step
 * updates each of its layers taking one of the following actions:
 * 1) no action (the layer remains unchanged - used to set different evolution speeds)
 * 2) internal evolution (the layer evolves according to some internal dynamics, defined by an EvolutionModel)
 * 3) external evolution (the layer imports nodes and edges from another layer)
 */

#ifndef UU_GENERATION_EVOLVE_H_
#define UU_GENERATION_EVOLVE_H_

#include <unordered_map>
#include "core/utils/random.hpp"
#include "core/exceptions/ElementNotFoundException.hpp"
#include "networks/_impl/containers/GenericObjectList.hpp"
#include "objects/Vertex.hpp"
#include "generation/EvolutionModel.hpp"

namespace uu {
namespace net {

//typedef int evolution_strategy;

//const int EVOLUTION_DEGREE=0;

/**
 * @brief Grows the input multilayer network.
 * @param net MLNetwork to grow
 * @param num_of_steps number of evolution steps
 * @param pr_no_event[] for each layer, the probability that during an evolution step the layer does not change
 * @param pr_internal_event[] for each layer, the probability that if something happens it is an internal evolution according to the evolution_model[] parameter
 * @param dependency[][] The (i,j) element of this matrix indicates the probability that, given an external evolution event, layer i will consider layer j as a potential candidate to import edges from
 * @param evolution_model[] for each layer, a specification of how the layer should evolve when an internal event happens
 **/
template <typename M>
void
evolve(
    M* net,
    const std::vector<std::string>& layer_names,
    const std::vector<double>& pr_internal_event,
    const std::vector<double>& pr_external_event,
    const std::vector<std::vector<double>>& dependency,
    const std::vector<EvolutionModel<M>*>& evolution_model,
    long num_of_steps
);



template <typename M>
void
evolve(
    M* net,
    const std::vector<std::string>& layer_names,
    const std::vector<double>& pr_internal_event,
    const std::vector<double>& pr_external_event,
    const std::vector<std::vector<double>>& dependency,
    const std::vector<EvolutionModel<M>*>& evolution_model,
    long num_of_steps
)
{

    std::unordered_map<std::string, size_t> layer_idx;

    for (size_t n=0; n<net->layers()->size(); n++)
    {
        layer_idx[net->layers()->at(n)->name] = n;
    }

    for (auto layer_name: layer_names)
    {
        if (layer_idx.find(layer_name) == layer_idx.end())
        {
            throw core::ElementNotFoundException("Layer " + layer_name + " not found");
        }
    }

    std::vector<double> pr_no_event;

    for (size_t i=0; i<pr_internal_event.size(); i++)
    {
        pr_no_event.push_back(1-pr_internal_event.at(i)-pr_external_event.at(i));
    }


    // Initialization
    std::vector<GenericObjectList<Vertex>> available_actors(net->layers()->size());

    for (size_t idx=0; idx<net->layers()->size(); idx++)
    {
        size_t n = layer_idx.at(layer_names.at(idx));
        available_actors[n] = GenericObjectList<Vertex>();

        for (auto actor: *net->actors())
        {
            available_actors[n].add(actor);
        }

        evolution_model[n]->init_step(net,net->layers()->at(n),available_actors[n]);

    }


    // Evolution
    for (long i=0; i<num_of_steps; i++)
    {
        for (size_t idx=0; idx<net->layers()->size(); idx++)
        {

            size_t n = layer_idx.at(layer_names.at(idx));

            auto target_layer = net->layers()->at(n);

            double dice = core::drand();

            if (dice < pr_no_event[n])
            {
                // DO NOTHING;
            }
            else if (dice < pr_internal_event[n]+pr_no_event[n] || pr_external_event[n]==0)
            {
                // INTERNAL EVOLUTION
                evolution_model[n]->internal_evolution_step(net,target_layer,available_actors[n]);

            }

            else
            {
                // EXTERNAL EVOLUTION
                // choose a layer from which to import edges.
                size_t layer_index = core::test(dependency[n]);
                auto external_layer = net->layers()->get(layer_names.at(layer_index));

                evolution_model[n]->external_evolution_step(net,target_layer,available_actors[n],external_layer);
            }
        }
    }
}

}
}

#endif
