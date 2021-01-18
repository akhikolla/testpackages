#ifndef UU_LAYOUT_CIRCULAR_H_
#define UU_LAYOUT_CIRCULAR_H_

#include <map>
#include "layout/XYZCoordinates.hpp"

namespace uu {
namespace net {

template <typename M>
std::map<std::pair<const Vertex*, const typename M::layer_type*>,XYZCoordinates>
circular(
    const M* mnet,
    double radius
);

}
}

#include "circular.ipp"

#endif
