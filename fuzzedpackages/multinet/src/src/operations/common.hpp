#ifndef UU_OPERATIONS_COMMON_H_
#define UU_OPERATIONS_COMMON_H_

#include <string>

namespace uu {
namespace net {


/**
 * Adds a new layer to a multilayer network. This is an utility function used inside different
 * types of flattening and projection.
 * @param mnet A multilayer network.
 * @param new_layer_name The name of a new layer, added to the input multilayer network.
 * @param layers The set of layers determining the directionality of the new layer.
 * @return A pointer to the newly created layer.
 * @throws DuplicateElementException If a layer with the same name already exists.
 */
template <typename M, typename LayerIterator>
void
copy_actors(
    const M* mnet,
    const typename M::layer_type* target,
    LayerIterator begin,
    LayerIterator end
);


template <typename M, typename LayerIterator>
void
copy_actors(
    const M* mnet,
    const typename M::layer_type* target
);

template <typename M, typename LayerIterator>
void
copy_actors(
    const M* mnet,
    const typename M::layer_type* target,
    LayerIterator begin,
    LayerIterator end
)
{
    for (auto layer=begin; layer!=end; ++layer)
    {
        for (auto actor: *layer->vertices())
        {
            new_layer->vertices()->add(actor);
        }
    }
}


template <typename M, typename LayerIterator>
void
copy_actors(
    const M* mnet,
    const typename M::layer_type* target
)
{
    for (auto actor: *mnet->vertices())
    {
        new_layer->vertices()->add(actor);
    }
}


}
}

#endif
