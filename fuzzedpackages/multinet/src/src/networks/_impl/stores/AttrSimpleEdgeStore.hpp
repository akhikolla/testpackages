/**
 * History:
 * - 2019.08.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_STORES_ATTRSIMPLEEDGESTORE_H_
#define UU_NET_STORES_ATTRSIMPLEEDGESTORE_H_

#include <memory>
#include "networks/_impl/stores/SimpleEdgeStore.hpp"
#include "core/stores/AttributeStore.hpp"

namespace uu {
namespace net {

/**
 * A SimpleEdgeStore allows to create, store, retrieve and erase a set of vertices.
 */
class
    AttrSimpleEdgeStore
    : public SimpleEdgeStore
{

  public:

    AttrSimpleEdgeStore(
        EdgeDir dir
    );


    core::AttributeStore<Edge>*
    attr(
    );


    const core::AttributeStore<Edge>*
    attr(
    ) const;

  private:

    std::unique_ptr<core::AttributeStore<Edge>> attributes_;

};


}
}

#endif
