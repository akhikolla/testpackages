#ifndef UU_OLAP_SEL_ENTRYITERATOR_H_
#define UU_OLAP_SEL_ENTRYITERATOR_H_

#include "core/exceptions/OutOfBoundsException.hpp"
#include "core/olap/selection/Indexes.hpp"
#include "core/olap/selection/IndexIterator.hpp"
#include <vector>

namespace uu {
namespace core {
namespace sel {

// @todo add const version

template <typename CUBE>
class EntryIterator
{
  private:

    CUBE* c_;

    IndexIterator idx_;

  public:

    /**
     * Creates a selector.
     */
    EntryIterator(
        CUBE* c,
        const std::vector<std::vector<size_t>>& indexes
    );

    /**
     * Creates a selector.
     */
    EntryIterator(
        CUBE* c,
        std::vector<Indexes>& indexes
    );

    /** Iterator over the objects in this collection. */
    class
        iterator
    {

      public:

        typedef ptrdiff_t difference_type; //almost always ptrdiff_t
        typedef typename CUBE::entry_type value_type; //almost always T
        typedef const typename CUBE::entry_type& reference; //almost always T& or const T&
        typedef typename CUBE::entry_type* pointer; //almost always T* or const T*
        typedef std::forward_iterator_tag iterator_category;

        iterator(
            CUBE* c,
            const IndexIterator::iterator& current
        );

        /*
        iterator(
            typename CUBE::container_type*
        );*/

        /** Return the object pointed by this iterator */
        typename CUBE::entry_type
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
            const EntryIterator<CUBE>::iterator& rhs
        );

        /** Checks if this iterator differs from the input one */
        bool
        operator!=(
            const EntryIterator<CUBE>::iterator& rhs
        );

      private:

        CUBE* c_;

        /** Entry currently pointed to by this iterator */
        IndexIterator::iterator current_;

    };

    /** Returns an iterator to the first object in the collection */
    EntryIterator<CUBE>::iterator
    begin(
    ) const;

    /** Returns an iterator after the last object in the collection */
    EntryIterator<CUBE>::iterator
    end(
    ) const;

};

/* TEMPLATE CODE */

template <typename CUBE>
EntryIterator<CUBE>::
EntryIterator(
    CUBE* c,
    const std::vector<std::vector<size_t>>& indexes
)
    : c_(c), idx_(IndexIterator(indexes))
{
    // @todo check null

}


template <typename CUBE>
EntryIterator<CUBE>::
EntryIterator(
    CUBE* c,
    std::vector<Indexes>& indexes
)
{
    // @todo check null
    c_ = c;

    size_t num_dimensions = c->dim().size();

    if (num_dimensions != indexes.size())
    {
        throw core::OutOfBoundsException("cube dimensions and indexes not matching in size");
    }

    std::vector<std::vector<size_t>> indexes2;

    auto size = c->size();

    for (size_t i = 0; i < num_dimensions; i++)
    {
        indexes[i].eval(size[i]);
        indexes2.push_back(std::vector<size_t>());

        while (indexes[i].has_next())
        {
            indexes2[i].push_back(indexes[i].next());
        }
    }

    idx_ = IndexIterator(indexes2);

}


template <typename CUBE>
typename EntryIterator<CUBE>::iterator
EntryIterator<CUBE>::
begin(
) const
{
    return iterator(c_, idx_.begin());
}


template <typename CUBE>
typename EntryIterator<CUBE>::iterator
EntryIterator<CUBE>::
end(
) const
{
    return iterator(c_, idx_.end());
}


template <typename CUBE>
EntryIterator<CUBE>::iterator::
iterator(
    CUBE* c,
    const IndexIterator::iterator& current
) :
    c_(c), current_(current)
{
}

/* @todo ??? check
template <typename CUBE>
EntryIterator<CUBE>::iterator::
iterator(
    typename CUBE::container_type* iter
) : c_(nullptr),
    indexes_({{},{}}),
current(iter)
{
}*/


template <typename CUBE>
typename EntryIterator<CUBE>::iterator
EntryIterator<CUBE>::iterator::
operator++(
)
{
    // PREFIX
    ++current_;

    return *this;
}


template <typename CUBE>
typename EntryIterator<CUBE>::iterator
EntryIterator<CUBE>::iterator::
operator++(
    int
)
{
    // POSTFIX
    auto tmp = *this;

    ++(*this);

    return tmp;
}


template <typename CUBE>
typename CUBE::entry_type
EntryIterator<CUBE>::iterator::
operator*(
)
{
    return c_->at(*current_);
}


template <typename CUBE>
bool
EntryIterator<CUBE>::iterator::
operator==(
    const EntryIterator<CUBE>::iterator& rhs
)
{
    return current_ == rhs.current_;
}


template <typename CUBE>
bool
EntryIterator<CUBE>::iterator::
operator!=(
    const EntryIterator<CUBE>::iterator& rhs
)
{
    return current_ != rhs.current_;
}

}
}
}

#endif
