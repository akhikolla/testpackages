#ifndef UU_OPERATIONS_FLATTEN_H_
#define UU_OPERATIONS_FLATTEN_H_

namespace uu {
namespace net {


/**
 * The "weighted" flattening approach adds an edge between a1 and a2 in the flattened network,
 * stored as a new layer, if the edge is present in any on the input layers, or combinations of them.
 * The weight on the edge indicates the number of layers where the edge is present.
 * Only intra-layer edges are considered.
 * @param mnet The multilayer network containing the layers to be merged.
 * @param new_layer_name The name of a new layer, added to the input multilayer network and obtained as a flattening of some of its layers.
 * @param layers The set of layers to be flattened.
 * @param force_directed If true, the flattened layer will contain directed edges. If false, it will contain directed edges only if at least one of the flattened layers do.
 * @param force_actors If true, all the actors in the multilayer network will be included in the flattened layer, even if they do not appear in any of the input layers.
 * @throws DuplicateElementException If a layer with the same name already exists.
 */
template <typename LayerIterator, typename W>
void
flatten_weighted(
    LayerIterator begin,
    LayerIterator end,
    W* target
);

/**
 * The "disjunctive" flattening approach, also known as or-flattening and unweighted-fattening,
 * adds an edge between a1 and a2 in the flattened network, stored as a new layer, if the edge
 * is present in any on the input layers, or combinations of them. There is no difference
 * in the result if an edge is present one or more times in the input layers.
 * Only intra-layer edges are considered.
 * @param mnet The multilayer network containing the layers to be merged.
 * @param new_layer_name The name of a new layer, added to the input multilayer network and obtained as a flattening of some of its layers.
 * @param layers The set of layers to be flattened.
 * @param force_directed If true, the flattened layer will contain directed edges. If false, it will contain directed edges only if at least one of the flattened layers do.
 * @param force_actors If true, all the actors in the multilayer network will be included in the flattened layer, even if they do not appear in any of the input layers.
 * @throws DuplicateElementException If a layer with the same name already exists.
 */
template <typename LayerIterator, typename G>
void
flatten_unweighted(
    LayerIterator begin,
    LayerIterator end,
    G* target
);

/* Temporary function for ASONAM paper :)
LayerSharedPtr
 flatten_obaida(
                              MLNetworkSharedPtr& mnet,
                              const std::string& new_layer_name,
                              const std::unordered_set<LayerSharedPtr>& layers,
                              bool force_directed,
                              bool force_actors,
                                double threshold);
 */

}
}

#include "flatten.ipp"

#endif
