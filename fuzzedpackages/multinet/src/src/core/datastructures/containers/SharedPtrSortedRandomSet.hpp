/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_CORE_DATASTRUCTURES_CONTAINERS_SHAREDPTRSORTEDRANDOMSET_H_
#define UU_CORE_DATASTRUCTURES_CONTAINERS_SHAREDPTRSORTEDRANDOMSET_H_

#include <memory>
#include "core/datastructures/containers/PtrSortedRandomSet.hpp"
#include "core/datastructures/containers/SharedPtrEQ.hpp"
#include "core/datastructures/containers/SharedPtrLT.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace core {

/**
 * A container to store shared_ptr's to objects.
 *
 * The typical use case is to store objects that can also be stored (and belong to) other
 * containers. Objects are deleted when they are no longer present in any container.
 */
template <typename E>
class SharedPtrSortedRandomSet
    : public PtrSortedRandomSet<E,std::shared_ptr<E>,SharedPtrLT<E>,SharedPtrEQ<E>>
{
  private:
    typedef PtrSortedRandomSet<E,std::shared_ptr<E>,SharedPtrLT<E>,SharedPtrEQ<E>> super;

  public:

    SharedPtrSortedRandomSet(
    ) : super() {};

    /**
     * Creates a sorted set optimized to store a pre-defined number of entries.
     * @param start_capacity the initial capacity for which the sorted set is optimized
     */
    SharedPtrSortedRandomSet(
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
        std::shared_ptr<E> element
    );

    /**
     * Adds an existing element to the store.
     *
     * @param element a pointer to an existing element
     * @return a pointer to the element, or nullptr if an element with the same key already exists in the store.
     **/
    virtual
    bool
    add(
        E * const element
    );

};



template <class E>
E*
SharedPtrSortedRandomSet<E>::
add(
    std::shared_ptr<E> element
)
{

    auto ptr = element.get();

    assert_not_null(ptr, "SharedPtrSortedRandomSet::add", "element");

    if (super::set.add(element))
    {
        return ptr;
    }

    else
    {
        return nullptr;
    }
}


template <class E>
bool
SharedPtrSortedRandomSet<E>::
add(
    E * const element
)
{

    assert_not_null(element, "SharedPtrSortedRandomSet::add", "element");

    return add(element->shared_from_this());
}


}
}

#endif

