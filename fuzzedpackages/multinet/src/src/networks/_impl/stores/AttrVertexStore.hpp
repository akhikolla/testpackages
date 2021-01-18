/**
 * History:
 * - 2019.08.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_STORES_ATTRVERTEXSTORE_H_
#define UU_NET_STORES_ATTRVERTEXSTORE_H_

#include <memory>
#include "networks/_impl/stores/VertexStore.hpp"
#include "core/stores/AttributeStore.hpp"

namespace uu {
namespace net {

/**
 * A VertexStore allows to create, store, retrieve and erase a set of vertices.
 */
class
    AttrVertexStore
    : public VertexStore
{

  public:

    AttrVertexStore();


    core::AttributeStore<Vertex>*
    attr(
    );


    const core::AttributeStore<Vertex>*
    attr(
    ) const;

  private:

    std::unique_ptr<core::AttributeStore<Vertex>> attributes_;

};


}
}

#endif
