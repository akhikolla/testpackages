/**
 * History:
 * - 2019.07.21 File created
 */

#ifndef UU_NET_DATASTRUCTURES_OLAP_SLICE_H_
#define UU_NET_DATASTRUCTURES_OLAP_SLICE_H_

#include <memory>
#include <string>
#include "networks/_impl/olap/VCube.hpp"
#include "core/olap/selection/EntryIterator.hpp"

namespace uu {
namespace net {

template <typename C>
core::sel::EntryIterator<C>
islice(
    C* cube,
    const std::vector<std::vector<size_t>>& indexes
)
{
    return core::sel::EntryIterator<C>(cube, indexes);
}


std::unique_ptr<VCube>
vslice(
    VCube* cube,
    const std::vector<std::vector<size_t>>& indexes,
    const std::string& name
);


}
}

#endif
