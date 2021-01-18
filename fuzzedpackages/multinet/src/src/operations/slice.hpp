
#ifndef UU_TNET_TRANSFORMATION_SLICE_H_
#define UU_TNET_TRANSFORMATION_SLICE_H_

#include <memory>
#include "networks/TemporalNetwork.hpp"
#include "networks/OrderedMultiplexNetwork.hpp"

namespace uu {
namespace net {


/**
* @brief Transforms a temporal network into an ordered multiplex network.
* @param tnet pointer to a temporal network
* @param num_partitions number of time slices
* @return a pointer to an ordered multiplex network
**/
std::unique_ptr<OrderedMultiplexNetwork>
slice_equal_time(
    const TemporalNetwork* tnet,
    size_t num_partitions
);

}
}

#endif
