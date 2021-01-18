/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_IO_MULTILAYERMETADATA_H_
#define UU_NET_IO_MULTILAYERMETADATA_H_

#include <vector>
#include <map>
#include <unordered_map>
#include "core/attributes/Attribute.hpp"
#include "networks/_impl/GraphType.hpp"

namespace uu {
namespace net {

struct MultilayerMetadata
{
    GraphType features;

    bool is_multiplex = true;

    std::unordered_map<std::string, GraphType> layers;

    std::map<std::pair<std::string,std::string>, bool> interlayer_dir;

    std::vector<core::Attribute> vertex_attributes;
    std::vector<core::Attribute> interlayer_edge_attributes;
    std::unordered_map<std::string, std::vector<core::Attribute>> intralayer_vertex_attributes;
    std::unordered_map<std::string, std::vector<core::Attribute>> intralayer_edge_attributes;
};

}
}

#endif
