/**
 * History:
 * - 2019.07.21 File created
 */

#include "networks/_impl/olap/slice.hpp"

namespace uu {
namespace net {

std::unique_ptr<VCube>
vslice(
    VCube* cube,
    const std::vector<std::vector<size_t>>& indexes,
    const std::string& name
)
{


    // @todo check not null

    // @todo check bounds

    // get dimensions and members from input cube
    auto res = cube->vslice(name, indexes);
    return res;

}
}
}

