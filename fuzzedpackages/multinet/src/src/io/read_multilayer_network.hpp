#ifndef UU_IO_READMULTILAYERNETWORK_H_
#define UU_IO_READMULTILAYERNETWORK_H_

#include <string>
#include <memory>
#include "networks/MultilayerNetwork.hpp"
#include "io/_impl/read_ml_common.hpp"

namespace uu {
namespace net {

std::unique_ptr<MultilayerNetwork>
read_attributed_homogeneous_multilayer_network(
    const std::string& infile,
    const std::string& name,
    char separator,
    bool align = false
);


template <>
Network*
read_layer(
    MultilayerNetwork* ml,
    const std::vector<std::string>& fields,
    size_t from_idx,
    size_t line_number
);

template <>
void
read_vertex(
    MultilayerNetwork* ml,
    const std::vector<std::string>& fields,
    const MultilayerMetadata& meta,
    size_t line_number
);

template <>
void
read_intralayer_vertex(
    MultilayerNetwork* ml,
    const std::vector<std::string>& fields,
    const MultilayerMetadata& meta,
    size_t line_number
);

template <>
void
read_intralayer_edge(
    MultilayerNetwork* ml,
    const std::vector<std::string>& fields,
    const MultilayerMetadata& meta,
    size_t line_number
);


template <>
void
read_interlayer_edge(
    MultilayerNetwork* ml,
    const std::vector<std::string>& fields,
    const MultilayerMetadata& meta,
    size_t line_number
);

}
}

#endif
