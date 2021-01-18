#ifndef UU_ALGORITHMS_BFS_H_
#define UU_ALGORITHMS_BFS_H_

#include <unordered_set>
#include <queue>
#include "core/exceptions/assert_not_null.hpp"
#include "objects/Vertex.hpp"
#include "objects/EdgeMode.hpp"

namespace uu {
namespace net {

/**
 * Performs a breadth-first search over a graph.
 */
template<typename G>
class BFS
{
  public:

    BFS(
        const G* g,
        const Vertex* v,
        EdgeMode mode
    );

    const Vertex*
    get_next();

  private:

    const G* g;
    EdgeMode mode;
    std::queue<const Vertex*> queue;
    std::unordered_set<const Vertex*> processed;

};

}
}

#import "BFS.ipp"

#endif
