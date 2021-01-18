/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_CORE_DATASTRUCTURES_CONTAINERS_UNIQUEPTRSORTEDRANDOMSET_H_
#define UU_CORE_DATASTRUCTURES_CONTAINERS_UNIQUEPTRSORTEDRANDOMSET_H_

#include <memory>
#include "core/datastructures/containers/PtrSortedRandomSet.hpp"
#include "core/datastructures/containers/UniquePtrEQ.hpp"
#include "core/datastructures/containers/UniquePtrLT.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace core {


/**
 * A container to store unique_ptr's to objects.
 *
 * The container is responsible for ownership: when the container is garbaged, all the
 * objects stored in the container are also deleted.
 */
template <typename E>
class
    UniquePtrSortedRandomSet :
    public PtrSortedRandomSet<E,std::unique_ptr<E>,UniquePtrLT<E>,UniquePtrEQ<E>>
{
  private:
    typedef PtrSortedRandomSet<E,std::unique_ptr<E>,UniquePtrLT<E>,UniquePtrEQ<E>> super;

  public:

    UniquePtrSortedRandomSet(
    ) : super() {};

    /**
     * Creates a sorted set optimized to store a pre-defined number of entries.
     * @param start_capacity the initial capacity for which the sorted set is optimized
     */
    UniquePtrSortedRandomSet(
        size_t start_capacity
    ) : super(start_capacity) {};


    /**
     * Inserts a new object in the collection.
     * @return true if KEY was not already present, false otherwise
     * (in which case the object is updated with the new value)
     */

    virtual
    E*
    add(
        std::unique_ptr<E> element
    );
};



template <class E>
E*
UniquePtrSortedRandomSet<E>::
add(
    std::unique_ptr<E> element
)
{
    auto res = element.get();

    assert_not_null(res, "UniquePtrSortedRandomSet::add", "element");

    if (super::set.add(std::move(element)))
    {
        return res;
    }

    else
    {
        return nullptr;
    }
}

}
}

#endif

