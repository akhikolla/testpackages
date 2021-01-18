#ifndef UU_OBJECTS_VERTEX_H_
#define UU_OBJECTS_VERTEX_H_

#include <string>
#include <memory>
#include <iostream>
#include "core/datastructures/objects/NamedObject.hpp"

namespace uu {
namespace net {

/**
 * A vertex in a graph.
 */
class
    Vertex :
    public core::NamedObject,
    public std::enable_shared_from_this<Vertex>
{

  public:

    typedef std::string key_type;

    /** Constructor. */
    Vertex(
        const std::string& name
    );

    /*static
    std::shared_ptr<const Vertex>
    create(
        const key_type& name
    );*/

    /** Output function, presenting a complete description of the vertex. */
    std::string
    to_string(
    ) const;

    const key_type key;

};

std::ostream&
operator<<(std::ostream& os, const Vertex& v);

}
}

#endif
