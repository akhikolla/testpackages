/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_MNET_COMMUNITY_VERTEXLAYERCOMMUNITY_H_
#define UU_MNET_COMMUNITY_VERTEXLAYERCOMMUNITY_H_

#include <utility>
#include "objects/Vertex.hpp"
#include "community/Community.hpp"

namespace uu {
namespace net {


template <typename G>
class
    VertexLayerCommunity :
    public Community<std::pair<const Vertex*, const G*>>
{};

}
}

#endif
