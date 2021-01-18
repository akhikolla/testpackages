/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NETCORE_DATASTRUCTURES_CONTAINERS_GENERICOBJECTLIST_H_
#define UU_NETCORE_DATASTRUCTURES_CONTAINERS_GENERICOBJECTLIST_H_


#include "core/datastructures/containers/SortedRandomSet.hpp"
#include <memory>

namespace uu {
namespace net {

template<typename O>
class
    GenericObjectList :
    public core::SortedRandomSet<const O*>
{
  public:
    static const std::unique_ptr<GenericObjectList<O>> empty;
};

template<typename O>
const std::unique_ptr<GenericObjectList<O>>
        GenericObjectList<O>::
        empty = std::make_unique<GenericObjectList<O>>();

}
}

#endif
