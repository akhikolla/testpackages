#ifndef UU_MEASURES_STRENGTH_H_
#define UU_MEASURES_STRENGTH_H_

#include <vector>
#include <algorithm>
#include "core/exceptions/assert_not_null.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include "core/utils/Counter.hpp"
#include "objects/Vertex.hpp"
#include "objects/EdgeMode.hpp"

namespace uu {
namespace net {

/**
 * Returns the maximum degree (\Delta) of a graph.
 * @param g input graph
 * @param mode to select IN, OUT, or INOUT degree
 * @return the degree of the vertex with the maximum degree in g
 */
template<typename G>
size_t
maximum_strength(
    const G* g,
    const EdgeMode mode = EdgeMode::INOUT
);


/**
 * Returns the minimum degree (\delta) of a graph.
 * @param g input graph
 * @param mode to select IN, OUT, or INOUT degree
 * @return the degree of the vertex with the minimum degree in g
 */
template<typename G>
size_t
minimum_strength(
    const G* g,
    const EdgeMode mode = EdgeMode::INOUT
);

/**
 * Returns the degree sequence of a graph.
 * @param g input graph
 * @param mode to select IN, OUT, or INOUT degree
 * @return a vector of length n with the ordered sequence of vertex degrees in g
 */
template<typename G>
std::vector<double>
strength_sequence(
    const G* g,
    const EdgeMode mode = EdgeMode::INOUT
);

/**
 * Returns the degree distribution of a graph.
 * @param g input graph
 * @param mode to select IN, OUT, or INOUT degree
 * @return a vector dd, where dd[i] is the number of vertices having degree i
 *
template<typename G>
std::vector<size_t>
degree_distribution(
    const G* g,
    const EdgeMode mode
);*/

/**
 * Returns the degree of a vertex.
 * @param g input graph
 * @param v input vertex
 * @param mode to select IN, OUT, or INOUT degree
 * @return the (mode-)degree of v in g
 */
template<typename G>
double
strength(
    const G* g,
    const Vertex* v,
    const EdgeMode mode = EdgeMode::INOUT
);


/** DEFINITIONS */


template<typename G>
double
strength(
    const G* g,
    const Vertex* v,
    const EdgeMode mode
)
{
    core::assert_not_null(g, "degree", "g");
    core::assert_not_null(g, "degree", "v");

    if (!g->is_weighted())
    {
        throw core::WrongParameterException("strength can only be computed on weighted graphs");
    }

    double s = 0;

    for (auto edge: *g->edges()->incident(v, mode))
    {
        auto w = g->get_weight(edge);

        if (!w.null)
        {
            s += w.value;

            if (edge->v1 == edge->v2)
            {
                // loops are counted twice
                s += w.value;
            }
        }
    }

    return s;
}


}
}

#endif
