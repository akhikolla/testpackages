#include "generation/empty_copy.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include "core/utils/Counter.hpp"

namespace uu {
namespace net {


template<typename G, typename ForwardIterator>
std::unique_ptr<G>
vertex_induced_subgraph(
    const G* g,
    ForwardIterator first,
    ForwardIterator last
)
{
    auto g_sub = empty_copy(g);

    core::Counter<const Edge*> edges;

    for (auto v_iter = first; v_iter != last; ++v_iter)
    {
        auto v = *v_iter;

        // check if the vertex belongs to g
        if (! g->vertices()->contains(v))
        {
            std::string err = "vertex " + v->to_string() + " not present in the input graph";
            throw core::WrongParameterException(err);
        }

        g_sub->vertices()->add(v);

        for (auto e: *g->edges()->incident(v))
        {
            edges.inc(e);
        }
    }

    for (auto pair: edges.map())
    {
        if (pair.second == 2)
        {
            g_sub->edges()->add(pair.first);
        }
    }

    return g_sub;
}


/**
 * Returns the subgraph induced by a set of edges.
 * @param g the input graph
 * @param first, last forward iterators to the initial and final
 * positions of the sequence of const Edge*'s.
 * The range used is [first,last).
 */
template<typename G, typename ForwardIterator>
std::unique_ptr<G>
edge_induced_subgraph(
    const G* g,
    ForwardIterator first,
    ForwardIterator last
)
{
    auto g_sub = empty_copy(g);

    for (auto e_iter = first; e_iter != last; ++e_iter)
    {
        auto e = *e_iter;

        // check if the vertex belongs to g
        if (! g->edges()->contains(e))
        {
            std::string err = "edge " + e->to_string() + " not present in the input graph";
            throw core::WrongParameterException(err);
        }

        g_sub->vertices()->add(e->v1);
        g_sub->vertices()->add(e->v2);
        g_sub->edges()->add(e);

    }

    return g_sub;
}


}
}

