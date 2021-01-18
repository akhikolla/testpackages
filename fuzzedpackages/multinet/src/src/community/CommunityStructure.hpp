/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_COMMUNITY_COMMUNITYSTRUCTURE_H_
#define UU_NET_COMMUNITY_COMMUNITYSTRUCTURE_H_

#include "core/datastructures/containers/UniquePtrSortedRandomSet.hpp"

namespace uu {
namespace net {

/**
 * A set of communities.
 */
template <typename C>
class CommunityStructure :
    public core::UniquePtrSortedRandomSet<C>
{
};

}
}

#endif

