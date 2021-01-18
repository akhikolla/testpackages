
#include "core/exceptions/assert_not_null.hpp"
#include "core/exceptions/ElementNotFoundException.hpp"
#include "core/exceptions/DuplicateElementException.hpp"

namespace uu {
namespace net {

template<typename G>
const Vertex*
edge_subdivision(
    G* g,
    const Edge* e,
    const std::string& vertex_name
)
{
    core::assert_not_null(g, "edge_subdivision", "g");
    core::assert_not_null(e, "edge_subdivision", "e");

    if (!g->edges()->contains(e))
    {
        throw core::ElementNotFoundException("edge " + e->to_string());
    }

    auto new_vertex = g->vertices()->add(vertex_name);

    if (!new_vertex)
    {
        throw core::DuplicateElementException("vertex " + vertex_name);
    }

    g->edges()->add(e->v1, new_vertex);
    g->edges()->add(new_vertex, e->v2);
    g->edges()->erase(e);

    return new_vertex;
}


}
}

