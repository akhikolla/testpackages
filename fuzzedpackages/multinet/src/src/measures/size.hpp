#ifndef UU_MEASURES_SIZE_H_
#define UU_MEASURES_SIZE_H_

namespace uu {
namespace net {

/**
 * Returns the number of edges in the graph.
 */
template<typename G>
size_t
size(
    const G* g
);

}
}

#include "size.ipp"

#endif
