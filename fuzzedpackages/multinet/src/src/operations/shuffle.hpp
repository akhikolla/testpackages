
#ifndef UU_OPERATIONS_SHUFFLE_H_
#define UU_OPERATIONS_SHUFFLE_H_

#include <memory>
#include "networks/TemporalNetwork.hpp"
#include "networks/OrderedMultiplexNetwork.hpp"

namespace uu {
namespace net {


/**
* @brief Shuffle edges in each layer.
* @param net a multilayer network
* @param num number of shufflings
**/
void
shuffle(
    uu::net::OrderedMultiplexNetwork* net,
    size_t num
);

}
}

#endif
