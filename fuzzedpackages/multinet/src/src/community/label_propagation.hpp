#ifndef UU_COMMUNITY_LABELPROPAGATIONSINGLE_H_
#define UU_COMMUNITY_LABELPROPAGATIONSINGLE_H_

#include <chrono>
#include "community/CommunityStructure.hpp"
#include "community/Community.hpp"
#include "community/_impl/cutils.hpp"
#include "objects/EdgeMode.hpp"
#include "objects/Vertex.hpp"
#include "core/utils/Counter.hpp"
#include "core/utils/random.hpp"

namespace uu {
namespace net {

template <typename G>
std::unique_ptr<CommunityStructure<Community<const Vertex*>>>
label_propagation(
    const G* net
);

}
}

#include "label_propagation.ipp"

#endif
