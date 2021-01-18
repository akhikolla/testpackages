/**
 * History:
 * - 2019.07.21 File created
 */

#include "networks/_impl/olap/resize.hpp"

namespace uu {
namespace net {

void
resize(
    VCube* cube,
    const std::string& dim,
    const std::string& member
)
{
    cube->resize(dim, member);
}

}
}

