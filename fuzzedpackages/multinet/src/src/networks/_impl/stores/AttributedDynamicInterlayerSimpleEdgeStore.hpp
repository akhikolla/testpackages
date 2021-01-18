/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_DATASTRUCTURES_STORES_ATTRIBUTEDDYNAMICINTERLAYERSIMPLEEDGESTORE_H_
#define UU_NET_DATASTRUCTURES_STORES_ATTRIBUTEDDYNAMICINTERLAYERSIMPLEEDGESTORE_H_

#include "networks/_impl/stores/DynamicInterlayerSimpleEdgeStore.hpp"
#include "networks/_impl/stores/Attributed.hpp"
#include "networks/Network.hpp"
#include "networks/_impl/stores/SimpleEdgeStore.hpp"
#include "core/datastructures/observers/Subject.hpp"

namespace uu {
namespace net {

template <typename V, typename L, typename A>
class
    AttributedDynamicInterlayerSimpleEdgeStore :
    public Attributed<A>,
    public DynamicInterlayerSimpleEdgeStore<V,L>
{
    typedef DynamicInterlayerSimpleEdgeStore<V,L> super;

    using Attributed<A>::attributes_;
    using super::attach;

  public:

    AttributedDynamicInterlayerSimpleEdgeStore(
        std::unique_ptr<A> attr
    );

};

template<typename V, typename L, typename A>
AttributedDynamicInterlayerSimpleEdgeStore<V,L,A>::
AttributedDynamicInterlayerSimpleEdgeStore(
    std::unique_ptr<A> attr
) :
    Attributed<A>(std::move(attr))
{
    attach(this->attr());
}

}
}

#endif
