/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_COMMUNITY_VERTEXCOMMUNITY_H_
#define UU_NET_COMMUNITY_VERTEXCOMMUNITY_H_

#include "community/Community.hpp"
#include "objects/Vertex.hpp"

namespace uu {
namespace net {

class
    VertexCommunity :
    public Community<const Vertex*>
{};

}
}

#endif
