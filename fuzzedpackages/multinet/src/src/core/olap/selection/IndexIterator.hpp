/**
 */

#ifndef UU_OLAP_SEL_INDEXITERATOR_H_
#define UU_OLAP_SEL_INDEXITERATOR_H_

#include <vector>
#include <cstddef>
#include <stddef.h>

namespace uu {
namespace core {
namespace sel {

class IndexIterator
{
  public:

    /**
     * Creates indexes (number of dimensions).
     */
    IndexIterator(
        const std::vector<size_t>& indexes
    );

    /**
     * Creates indexes.
     */
    IndexIterator(
        const std::vector<std::vector<size_t>>& indexes
    );


    /** Iterator over the objects in this collection. */
    class
        iterator
    {

      public:

        typedef ptrdiff_t difference_type; //almost always ptrdiff_t
        typedef std::vector<size_t> value_type; //almost always T
        typedef const std::vector<size_t>& reference; //almost always T& or const T&
        typedef std::vector<size_t>* pointer; //almost always T* or const T*
        typedef std::forward_iterator_tag iterator_category;


        iterator(
            const std::vector<std::vector<size_t>>& indexes
        );

        // @todo check encapsulation
        iterator(
            const std::vector<std::vector<size_t>>& indexes,
            const std::vector<size_t>& current
        );

        /** Return the object pointed by this iterator */
        std::vector<size_t>
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
            const IndexIterator::iterator& rhs
        );

        /** Checks if this iterator differs from the input one */
        bool
        operator!=(
            const IndexIterator::iterator& rhs
        );

      private:

        /* selection expressed as indexes of each dimension */
        std::vector<std::vector<size_t>> indexes_;

        /** Entry currently pointed to by this iterator */
        std::vector<size_t> current_;

        /** Order in which to increment the indexes */
        std::vector<size_t> pivot_;
    };

    iterator
    begin(
    ) const;

    iterator
    end(
    ) const;

  private:

    /* selection expressed as indexes of each dimension */
    std::vector<std::vector<size_t>> indexes_;

};

}
}
}

#endif
