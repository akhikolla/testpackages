
#include "core/exceptions/assert_not_null.hpp"
#include "core/exceptions/ElementNotFoundException.hpp"
#include "core/exceptions/DuplicateElementException.hpp"

namespace uu {
namespace net {

template<typename G>
const Vertex*
edge_contraction(
    G* g,
    const Edge* e,
    const std::string& vertex_name
)
{

    core::assert_not_null(g, "edge_contraction", "g");
    core::assert_not_null(e, "edge_contraction", "e");

    if (!g->edges()->contains(e))
    {
        throw core::ElementNotFoundException("edge " + e->to_string());
    }

    auto new_vertex = g->vertices()->add(vertex_name);

    if (!new_vertex)
    {
        throw core::DuplicateElementException("vertex " + vertex_name);
    }

    if (!g->is_directed())
    {
        for (auto neigh_v1: *g->edges()->neighbors(e->v1))
        {
            if (neigh_v1 == e->v2)
            {
                continue;
            }

            g->edges()->add(neigh_v1, new_vertex);
        }

        for (auto neigh_v2: *g->edges()->neighbors(e->v2))
        {
            if (neigh_v2 == e->v1)
            {
                continue;
            }

            g->edges()->add(new_vertex, neigh_v2);
        }
    }

    else
    {
        for (auto neigh_v1: *g->edges()->neighbors(e->v1, EdgeMode::IN))
        {
            if (neigh_v1 == e->v2)
            {
                continue;
            }

            g->edges()->add(neigh_v1, new_vertex);
        }

        for (auto neigh_v1: *g->edges()->neighbors(e->v1, EdgeMode::OUT))
        {
            if (neigh_v1 == e->v2)
            {
                continue;
            }

            g->edges()->add(new_vertex, neigh_v1);
        }

        for (auto neigh_v2: *g->edges()->neighbors(e->v2, EdgeMode::IN))
        {
            if (neigh_v2 == e->v1)
            {
                continue;
            }

            g->edges()->add(neigh_v2, new_vertex);
        }

        for (auto neigh_v2: *g->edges()->neighbors(e->v2, EdgeMode::OUT))
        {
            if (neigh_v2 == e->v1)
            {
                continue;
            }

            g->edges()->add(new_vertex, neigh_v2);
        }
    }

    auto v1 = e->v1;
    auto v2 = e->v2;
    g->vertices()->erase(v1);
    g->vertices()->erase(v2);

    return new_vertex;
}


}
}

