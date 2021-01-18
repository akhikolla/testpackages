/**
 * This module defines a generic network co-evolution process.
 *
 * The function "evolve" takes a multilayer network as input and at every step
 * updates each of its layers taking one of the following actions:
 * 1) no action (the layer remains unchanged - used to set different evolution speeds)
 * 2) internal evolution (the layer evolves according to some internal dynamics, defined by an EvolutionModel)
 * 3) external evolution (the layer imports nodes and edges from another layer)
 */

#ifndef UU_GENERATION_VERTEXSELECTION_H_
#define UU_GENERATION_VERTEXSELECTION_H_


namespace uu {
namespace net {

template <typename G>
const Vertex*
choice_uniform(
    const G* g
);

template <typename G>
const Vertex*
choice_common_neighbors(
    const G* g1,
    const G* g2
);

template <typename G>
const Vertex*
choice_degree(
    const G* layer
);


}
}

#endif
