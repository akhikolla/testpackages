#ifndef UU_IO_WRITEMULTILAYERNETWORK_H_
#define UU_IO_WRITEMULTILAYERNETWORK_H_

#include <string>
#include <memory>
#include "networks/MultilayerNetwork.hpp"
#include "core/utils/string.hpp"
#include "io/_impl/read_common.hpp"

namespace uu {
namespace net {

template <typename LayerIterator>
void
write_attributed_homogeneous_multilayer_network(
    const MultilayerNetwork* net,
    LayerIterator begin,
    LayerIterator end,
    const std::string& outfile,
    char separator
);

template <typename LayerIterator>
void
write_graphml(
    const MultilayerNetwork* mnet,
    LayerIterator begin,
    LayerIterator end,
    const std::string& path,
    bool merge_actors,
    bool include_all_actors
);

}
}

#include "write_multilayer_network.ipp"

#endif
