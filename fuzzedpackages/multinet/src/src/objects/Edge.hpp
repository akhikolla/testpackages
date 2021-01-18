#ifndef UU_OBJECTS_EDGE_H_
#define UU_OBJECTS_EDGE_H_

#include <string>
#include <memory>
#include <iostream>
#include "objects/Vertex.hpp"
#include "objects/EdgeDir.hpp"
#include "core/datastructures/objects/Object.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace net {

/**
 * An edge between two vertices.
 */
class
    Edge :
    public core::Object,
    public std::enable_shared_from_this<Edge>
{

  public:

    /** Constructor. */
    Edge(
        const Vertex* v1,
        const Vertex* v2,
        EdgeDir dir
    );

    static
    std::shared_ptr<Edge>
    create(
        const Vertex* v1,
        const Vertex* v2,
        EdgeDir dir
    );

    /** Output function, presenting a complete description of the edge. */
    std::string
    to_string(
    ) const;

    /** The vertex at the first end of this edge. */
    const Vertex* v1;

    /** The vertex at the second end of this edge. */
    const Vertex* v2;

    /** Edge directionality. */
    const EdgeDir dir;

};

std::ostream&
operator<<(
    std::ostream& os,
    const Edge& e
);


}
}

#endif
