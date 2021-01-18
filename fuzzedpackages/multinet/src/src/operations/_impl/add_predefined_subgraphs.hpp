#ifndef UU_NET_OPERATIONS_ADDPREDEFINEDSUBGRAPHS_H_
#define UU_NET_OPERATIONS_ADDPREDEFINEDSUBGRAPHS_H_

namespace uu {
namespace net {


/**
 * Adds n new vertices to a graph.
 *
 * G must have a method vertices() returning a vertex store.
 */
template<typename G>
std::vector<const Vertex*>
add_vertices(
    G* g,
    size_t n,
    const std::string& base_vertex_name = "v"
);

/**
 * Adds n new vertices to a graph, if no vertices with the same names exist.
 * @param first, last forward iterators to the initial and final
 * positions of the sequence of const Vertex*'s.
 * The range used is [first,last).
 *
 * G must have a method vertices() returning a vertex store.
 */
template<typename G, typename ForwardIterator>
void
add_vertices(
    G* g,
    ForwardIterator first,
    ForwardIterator last
);

/**
 * Adds n new vertices, all adjacent to each other, if no vertices with the same names exist.
 * If the graph is directed, arcs in both directions are created for each pair of vertices.
 * No loops are created.
 */
template<typename G>
void
add_complete_subgraph(
    G* g,
    size_t n,
    const std::string& base_vertex_name = "v"
);

/**
 * Adds two partite sets V1 and V2 where all vertices in V1 are all adjacent
 * to all vertices in V2, if no vertices with the same names exist.
 * If the graph is directed, arcs in both directions are created for each pair of vertices.
 */
template<typename G>
void
add_complete_bipartite_subgraph(
    G* g,
    size_t n1,
    size_t n2,
    const std::string& base_vertex_name1 = "v",
    const std::string& base_vertex_name2 = "u"
);


/**
 * Adds n new vertices with edges forming a path, if no vertices with the same names exist.
 */
template<typename G>
void
add_path(
    G* g,
    size_t n,
    const std::string& base_vertex_name = "v"
);

/**
 * Adds n new vertices with edges forming a circle, if no vertices with the same names exist.
 */
template<typename G>
void
add_cycle(
    G* g,
    size_t n,
    const std::string& base_vertex_name = "v"
);


/**
 * Adds n new vertices made of a C_{n-1}
 * cycle graph plus a vertex adjacent to all the other n-1 vertices,
 * if no vertices with the same names exist.
 */
template<typename G>
void
add_wheel(
    G* g,
    size_t n,
    const std::string& base_vertex_name = "v"
);

}
}

#include "add_predefined_subgraphs.ipp"

#endif
