/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_IO_GRAPHMETADATA_H_
#define UU_NET_IO_GRAPHMETADATA_H_

#include <vector>
#include "core/attributes/Attribute.hpp"
#include "networks/_impl/GraphType.hpp"

namespace uu {
namespace net {

struct GraphMetadata
{
    GraphType features;
    std::vector<core::Attribute> vertex_attributes;
    std::vector<core::Attribute> edge_attributes;
};

}
}

#endif
