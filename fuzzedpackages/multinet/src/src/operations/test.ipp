#include <unordered_set>
#include <unordered_map>
#include <queue>
#include "core/exceptions/assert_not_null.hpp"
#include "objects/Vertex.hpp"
#include "objects/EdgeMode.hpp"

namespace uu {
namespace net {


template<typename G>
bool
is_bipartite(
    const G*  g
)
{

    core::assert_not_null(g, "is_bipartite", "g");

    // Each vertex can be assigned to one of two partitions: A (represented by value true)
    // or B (represented by value false)
    std::unordered_map<const Vertex*, bool> partition;

    // We perform a BFS from each vertex, using a queue

    std::queue<const Vertex*> queue;

    for (auto v: *g->vertices())
    {
        auto p = partition.find(v);

        // We process vertex v only if it has not been assigned to any partition yet

        if (p == partition.end())
        {
            queue.push(v);
            partition[v] = true;

            while (queue.size()>0)
            {
                auto u = queue.front();
                queue.pop();
                bool part_u = partition[u];

                // We try to assign all neighbors of u to the other partition (!part_u)

                for (auto n: *g->edges()->neighbors(u, EdgeMode::INOUT))
                {
                    auto p_n = partition.find(n);

                    if (p_n == partition.end())
                    {
                        // Vertex not assigned to any partition yet, so we assign it

                        partition[n] = !part_u;
                        queue.push(n);
                    }

                    else
                    {
                        // Vertex already in a partition. We check if it is the right one.
                        bool part_n = p_n->second;

                        if (part_n == part_u)
                        {
                            // Two neighbors in the same partition: test failed.

                            return false;
                        }
                    }
                }
            }
        }
    }

    return true;
}

}
}
