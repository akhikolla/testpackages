#ifndef UU_ALGORITHMS_DFS_H_
#define UU_ALGORITHMS_DFS_H_

#include <unordered_set>
#include <stack>
#include "core/exceptions/assert_not_null.hpp"
#include "objects/Vertex.hpp"
#include "objects/EdgeMode.hpp"

namespace uu {
namespace net {

/**
 * Performs a depth-first search over a graph.
 */
template<typename G>
class DFS
{
  public:

    DFS(
        const G* g,
        const Vertex* v,
        EdgeMode mode
    );

    const Vertex*
    get_next();

  private:

    const G* g;
    EdgeMode mode;
    std::stack<const Vertex*> stack;
    std::unordered_set<const Vertex*> processed;

};

}
}

#import "DFS.ipp"

#endif
