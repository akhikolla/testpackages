#ifndef UU_MEASURES_DEGREE_H_
#define UU_MEASURES_DEGREE_H_

#include <vector>
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
maximum_degree(
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
minimum_degree(
    const G* g,
    const EdgeMode mode = EdgeMode::INOUT
);


/**
 * Returns the average degree of a graph.
 * @param g input graph
 * @param mode to select IN, OUT, or INOUT degree
 * @return the average degree of vertices in g
 */
template<typename G>
double
average_degree(
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
std::vector<size_t>
degree_sequence(
    const G* g,
    const EdgeMode mode = EdgeMode::INOUT
);

/**
 * Returns the degree distribution of a graph.
 * @param g input graph
 * @param mode to select IN, OUT, or INOUT degree
 * @return a vector dd, where dd[i] is the number of vertices having degree i
 */
template<typename G>
std::vector<size_t>
degree_distribution(
    const G* g,
    const EdgeMode mode = EdgeMode::INOUT
);

/**
 * Returns the degree of a vertex.
 * @param g input graph
 * @param v input vertex
 * @param mode to select IN, OUT, or INOUT degree
 * @return the (mode-)degree of v in g
 */
template<typename G>
size_t
degree(
    const G* g,
    const Vertex* v,
    const EdgeMode mode = EdgeMode::INOUT
);

}
}

#include "degree.ipp"

#endif
