/**
 * This header defines an iterator producing object names.
 *
 * History:
 * - 2018.01.01 file imported from version 1.0 of the multinet library
 */

#ifndef UU_CORE_UTILS_NAMES_H_
#define UU_CORE_UTILS_NAMES_H_

#include <cstddef>
#include "stddef.h"
#include <string>

namespace uu {
namespace core {


class NameIterator
{
  public:

    /**
     * Creates names.
     */
    NameIterator(
        const std::string& prefix,
        size_t num_names
    );


    /** Iterator over the names. */
    class
        iterator
    {

      public:

        typedef ptrdiff_t difference_type;
        typedef std::string value_type;
        typedef const std::string& reference;
        typedef std::string* pointer;
        typedef std::forward_iterator_tag iterator_category;

        iterator(
            const std::string& prefix,
            size_t num_names,
            size_t num_digits
        );

        iterator(
            const std::string& prefix,
            size_t num_names,
            size_t num_digits,
            size_t current
        );

        /** Returns the name pointed by this iterator */
        std::string
        operator*(
        );

        /** Moves the iterator to the next name (prefix) */
        iterator
        operator++(
        );

        /** Moves the iterator to the next name */
        iterator
        operator++(
            int
        );

        /** Checks if this iterator equals the input one */
        bool
        operator==(
            const NameIterator::iterator& rhs
        );

        /** Checks if this iterator differs from the input one */
        bool
        operator!=(
            const NameIterator::iterator& rhs
        );

      private:

        std::string prefix_;

        size_t num_names_;

        size_t current_;

        size_t num_digits_;

    };

    iterator
    begin(
    ) const;

    iterator
    end(
    ) const;

  private:

    std::string prefix_;

    size_t num_names_;

    size_t num_digits_;
};

}
}

#endif
