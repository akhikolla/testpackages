
#include "generation/empty_copy.hpp"

namespace uu {
namespace net {

template<typename G>
std::unique_ptr<G>
graph_complement(
    const G* g,
    const std::string& name
)
{
    std::unique_ptr<G> res = empty_copy(g, name);

    for (auto vertex: *g->vertices())
    {
        res->vertices()->add(vertex);
    }

    for (auto v1: *g->vertices())
    {
        for (auto v2: *g->vertices())
        {
            if (!res->allows_loops() && (v1==v2))
            {
                continue;
            }

            if (!res->is_directed() && (v1>v2))
            {
                continue;
            }

            if (!g->edges()->get(v1, v2))
            {
                res->edges()->add(v1, v2);
            }
        }
    }

    return res;
}

}
}

