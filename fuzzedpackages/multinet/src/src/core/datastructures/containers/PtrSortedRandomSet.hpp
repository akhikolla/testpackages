
/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_CORE_DATASTRUCTURES_CONTAINERS_PTRSORTEDRANDOMSET_H_
#define UU_CORE_DATASTRUCTURES_CONTAINERS_PTRSORTEDRANDOMSET_H_

#include <memory>
#include "core/datastructures/containers/SortedRandomSet.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace core {

/**
 * A generic container to store pointers to objects.
 * This class provides common functionality used by UniquePtrSortedRandomSet
 * and SharedPtrSortedRandomSet
 *
 * This class requires to specify both the type of the pointed object (E) and
 * the pointer type (inclusive of the type of the pointed object).
 * A template template would be cleaner, but this construct is still not fully supported by all
 * compilers.
 */
template <typename E, class PTR, typename PtrLT, typename PtrEQ>
class
    PtrSortedRandomSet
{

  protected:
    SortedRandomSet<PTR> set;

  public:

    PtrSortedRandomSet(
    );

    /**
     * Creates a sorted set optimized to store a pre-defined number of entries.
     * @param start_capacity the initial capacity for which the sorted set is optimized
     */
    PtrSortedRandomSet(
        size_t start_capacity
    );


    /** Iterator over the objects in this collection. */
    class
        iterator
    {
        typedef std::forward_iterator_tag iterator_category;

      public:

        iterator();

        iterator(
            typename SortedRandomSet<PTR>::iterator
        );

        /** Return the object pointed by this iterator */
        E*
        operator*(
        );

        /** Moves the iterator to the next object in the collection (prefix) */
        iterator
        operator++(
        );

        /** Moves the iterator to the next object in the collection (postfix) */
        iterator
        operator++(
            int
        );

        /** Checks if this iterator equals the input one */
        bool
        operator==(
            const PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator& rhs
        );

        /** Checks if this iterator differs from the input one */
        bool
        operator!=(
            const PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator& rhs
        );

      protected:
        typename SortedRandomSet<PTR>::iterator it;
    };

    /** Returns an iterator to the first object in the collection */
    PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::
    iterator
    begin(
    ) const;

    /** Returns an iterator after the last object in the collection */
    PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::
    iterator
    end(
    ) const;

    /** Returns the number of objects in the collection */
    size_t
    size(
    ) const;

    /** Returns true if an object with the input id is present in the collection */
    bool
    contains(
        const E*
    ) const;

    /** Returns the position of the input value in the collection, or -1 */
    int
    index_of(
        const E*
    ) const;

    /** Returns the object at the given position in the collection.
     * @throw ElementNotFoundException if the index is outside the bounds on the set
     */
    E*
    at(
        size_t
    ) const;

    /** Returns a random object, uniform probability */
    E*
    get_at_random(
    ) const;

    /**
     * Inserts a new object in the collection.
     * @return true if KEY was not already present, false otherwise
     * (in which case the object is updated with the new value)
     */

    virtual
    E*
    add(
        PTR element
    ) = 0;

    /**
     * Erases an existing vertex from the store.
     *
     * @param vertex a pointer to the vertex to be deleted
     * @return true if the object has been erased, false if it was not present.
     **/
    virtual
    bool
    erase(
        E * const element
    );

};


template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::
PtrSortedRandomSet(
)
{
}

template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::
PtrSortedRandomSet(
    size_t start_capacity
) : set(start_capacity)
{
}


template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
bool
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::
erase(
    E * const element
)
{
    core::assert_not_null(element, "erase", "e");

    return set.template erase<E*,PtrLT,PtrEQ>(element);
}




template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
typename PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::
begin(
) const
{
    return iterator(set.begin());
}

template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
typename PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::
end(
) const
{
    return iterator(set.end());
}

template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
E*
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator::operator*(
)
{
    return (*it).get();
}


template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator::
iterator(
    typename SortedRandomSet<PTR>::iterator iter
) :
    it(iter)
{
}


template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
typename PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator::
operator++(
)
{
    // PREFIX
    ++it;
    return *this;
}

template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
typename PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator::
operator++(
    int
)
{
    // POSTFIX
    PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator tmp(it);
    ++it;
    return tmp;
}

template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
bool
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator::
operator==(
    const PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator& rhs
)
{
    return it == rhs.it;
}

template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
bool
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator::
operator!=(
    const PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::iterator& rhs
)
{
    return it != rhs.it;
}

template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
size_t
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::size(
) const
{
    return set.size();
}

template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
bool
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::
contains(
    const E* search_value
) const
{

    core::assert_not_null(search_value, "contains", "search_value");

    return set.template contains<const E*,PtrLT,PtrEQ>(search_value);
}

template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
int
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::
index_of(
    const E* search_value
) const
{

    core::assert_not_null(search_value, "get_index", "search_value");

    return set.template index_of<const E*,PtrLT,PtrEQ>(search_value);
}

template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
E*
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::
at(
    size_t pos
) const
{
    return set.at(pos).get();
}

template <typename E, typename PTR, typename PtrLT, typename PtrEQ>
E*
PtrSortedRandomSet<E, PTR, PtrLT, PtrEQ>::
get_at_random(
) const
{
    return set.get_at_random().get();
}


}
}

#endif

