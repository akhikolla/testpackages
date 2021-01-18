#ifndef UU_IO_READPROBABILISTICNETWORK_H_
#define UU_IO_READPROBABILISTICNETWORK_H_

#include <string>
#include <memory>
#include "networks/ProbabilisticNetwork.hpp"
#include "io/_impl/read_common.hpp"

namespace uu {
namespace net {

std::unique_ptr<ProbabilisticNetwork>
read_probabilistic_network(
    const std::string& infile,
    const std::string& name,
    char separator
);

template <>
void
read_vertex(
    ProbabilisticNetwork* g,
    const std::vector<std::string>& fields,
    const std::vector<core::Attribute>& vertex_attributes,
    size_t line_number
);


template <>
void
read_edge(
    ProbabilisticNetwork* g,
    const std::vector<std::string>& fields,
    const std::vector<core::Attribute>& edge_attributes,
    size_t line_number
);


}
}

#endif
