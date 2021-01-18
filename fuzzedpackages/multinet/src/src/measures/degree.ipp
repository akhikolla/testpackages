
#include <algorithm>
#include "core/exceptions/assert_not_null.hpp"
#include "core/utils/Counter.hpp"
#include "measures/order.hpp"
#include "measures/size.hpp"

namespace uu {
namespace net {

template<typename G>
size_t
maximum_degree(
    const G* g,
    const EdgeMode mode
)
{
    core::assert_not_null(g, "maximum_degree", "g");
    size_t max = 0;
    size_t d;

    for (auto v: *g->vertices())
    {
        d=degree(g, v, mode);

        if (d > max)
        {
            max = d;
        }
    }

    return max;
}





template<typename G>
size_t
minimum_degree(
    const G* g,
    const EdgeMode mode
)
{
    core::assert_not_null(g, "minimum_degree", "g");
    size_t min = 0;
    size_t d;
    bool first = true;

    for (auto v: *g->vertices())
    {
        d=degree(g, v, mode);

        if (first)
        {
            min = d;
            first = false;
        }

        else if (d < min)
        {
            min = d;
        }
    }

    return min;
}




template<typename G>
std::vector<size_t>
degree_sequence(
    const G* g,
    const EdgeMode mode
)
{
    core::assert_not_null(g, "degree_sequence", "g");

    size_t order = g->vertices()->size();
    std::vector<size_t> res;
    res.reserve(order);

    size_t d;

    for (auto v: *g->vertices())
    {
        d=degree(g, v, mode);
        res.push_back(d);
    }

    std::sort(res.begin(), res.end(), std::greater<size_t>());
    return res;
}




template<typename G>
std::vector<size_t>
degree_distribution(
    const G* g,
    const EdgeMode mode
)
{
    core::assert_not_null(g, "degree_distribution", "g");
    core::Counter<size_t> dd;

    size_t max = 0;
    size_t d;

    for (auto v: *g->vertices())
    {
        d = degree(g, v, mode);
        dd.inc(d);

        if (d > max)
        {
            max = d;
        }
    }

    std::vector<size_t> res;
    res.reserve(max+1);

    for (d = 0; d<=max; d++)
    {
        res.push_back(dd.count(d));
    }

    return res;
}


template<typename G>
double
average_degree(
    const G* g,
    const EdgeMode mode
)
{
    core::assert_not_null(g, "average_degree", "g");

    if (g->is_directed() && mode != EdgeMode::INOUT)
    {
        return size(g) / (double)order(g);
    }

    else
    {
        return 2.0 * size(g) / order(g);
    }
}

template<typename G>
size_t
degree(
    const G* g,
    const Vertex* v,
    const EdgeMode mode
)
{
    core::assert_not_null(g, "degree", "g");
    core::assert_not_null(g, "degree", "v");
    auto inc = g->edges()->incident(v, mode);
    auto d = inc->size();

    if (g->allows_loops())
    {
        for (auto e: *inc)
        {
            if (e->v1 == e->v2)
            {
                d++;    // for loops
            }
        }
    }

    return d;
}

}
}

