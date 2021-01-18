#ifndef UU_COMMUNITY_OPTIMIZESLICES_H_
#define UU_COMMUNITY_OPTIMIZESLICES_H_

#include "networks/TemporalNetwork.hpp"
#include <vector>

namespace uu {
namespace net {

std::vector<double>
optimize_slices(
    const TemporalNetwork* original_net,
    size_t max_slices);

}
}

#endif
