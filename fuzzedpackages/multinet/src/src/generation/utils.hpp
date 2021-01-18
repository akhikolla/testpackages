#ifndef UU_NET_GENERATION_UTILS_H_
#define UU_NET_GENERATION_UTILS_H_


#include <utility>
#include "objects/Vertex.hpp"

namespace uu {
namespace net {


template <typename G>
std::pair<const Vertex*, const Vertex*>
get_vertex_pair(
    G* g,
    bool allow_loops = false
);

}
}

#include "utils.ipp"

#endif
