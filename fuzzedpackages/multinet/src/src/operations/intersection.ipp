#include "core/exceptions/assert_not_null.hpp"
#include "core/exceptions/OperationNotSupportedException.hpp"
#include "generation/empty_copy.hpp"

namespace uu {
namespace net {

/**
 * Returns the intersection of two graphs.
 *
 * The operation is only allowed if both graphs are directed or both are undirected.
 *
 * Only vertices and edges are included in the new graph, not attributes.
 *
 * @param g1, g2 input graphs
 */
template<typename G>
std::unique_ptr<G>
graph_intersection(
    const G* g1,
    const G* g2,
    const std::string& name
)
{
    core::assert_not_null(g1, "graph_intersection", "g1");
    core::assert_not_null(g2, "graph_intersection", "g2");

    if (g1->is_directed() != g2->is_directed())
    {
        std::string err = "intersection between directed and undirected graphs";
        throw core::OperationNotSupportedException(err);
    }

    std::unique_ptr<G> res = empty_copy(g1, name);

    for (auto vertex: *g1->vertices())
    {
        if (g2->vertices()->contains(vertex))
        {
            res->vertices()->add(vertex);
        }
    }


    for (auto edge: *g1->edges())
    {
        if (g2->edges()->contains(edge))
        {
            res->edges()->add(edge);
        }
    }

    return res;
}

}
}

