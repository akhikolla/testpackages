#ifndef UU_CREATION_STANDARDGRAPHS_H_
#define UU_CREATION_STANDARDGRAPHS_H_

#include "networks/Network.hpp"

#include <memory>

namespace uu {
namespace net {


/**
 * Returns a network with n vertices and no edges.
 */
std::unique_ptr<Network>
null_graph(
    size_t n,
    EdgeDir dir = EdgeDir::UNDIRECTED,
    bool allows_loops = false
);


/**
 * Returns a network K_n with n vertices, all adjacent to each other.
 * If the graph is directed, arcs in both directions are created for each pair of vertices.
 * No loops are created.
 */
std::unique_ptr<Network>
complete_graph(
    size_t n,
    EdgeDir dir = EdgeDir::UNDIRECTED
);


/**
 * Returns a bipartite network K_{n1,n2} with two partite sets V and U
 * where all vertices in V are all adjacent to all vertices in U.
 * If the graph is directed, arcs in both directions are created for each pair of vertices.
 */
std::unique_ptr<Network>
complete_bipartite_graph(
    size_t n1,
    size_t n2,
    EdgeDir dir = EdgeDir::UNDIRECTED
);


/**
 * Returns a network P_n with n vertices with edges forming a path.
 */
std::unique_ptr<Network>
path_graph(
    size_t n,
    EdgeDir dir = EdgeDir::UNDIRECTED
);

/**
 * Returns a network C_n with n vertices with edges forming a circle.
 */
std::unique_ptr<Network>
cycle_graph(
    size_t n,
    EdgeDir dir = EdgeDir::UNDIRECTED
);


/**
 * Returns a network W_n with n vertices made of a C_{n-1}
 * cycle graph plus a vertex adjacent to all the other n-1 vertices.
 * The graph is undirected.
 */
std::unique_ptr<Network>
wheel_graph(
    size_t n
);

}
}


#endif
