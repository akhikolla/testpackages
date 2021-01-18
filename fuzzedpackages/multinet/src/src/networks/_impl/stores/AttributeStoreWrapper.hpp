/**
 * History:
 * - 2018.03.09 file created, with code taken from other existing files.
 */

#ifndef UU_NET_DATASTRUCTURES_STORES_ATTRIBUTESTOREWRAPPER_H_
#define UU_NET_DATASTRUCTURES_STORES_ATTRIBUTESTOREWRAPPER_H_

#include <memory>
#include "core/stores/AttributeStore.hpp"

namespace uu {
namespace net {

/**
 * A class to handle attributes to be associated to objects of type OT.
 *
 * OT can currently be Vertex or Edge.
 */
template <typename OT>
class AttributeStoreWrapper
{

  public:

    /**
     * Constructor.
     */
    AttributeStoreWrapper(std::unique_ptr<core::AttributeStore<OT>>);


  protected:

    std::unique_ptr<core::AttributeStore<OT>> attr_;

};

template <typename OT>
AttributeStoreWrapper<OT>::
AttributeStoreWrapper(
    std::unique_ptr<core::AttributeStore<OT>> attr
) : attr_(std::move(attr))
{
}



} // namespace net
} // namespace uu

#endif
