/**
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_EXCEPTIONS_DUPLICATEELEMENTEXCEPTION_H_
#define UU_CORE_EXCEPTIONS_DUPLICATEELEMENTEXCEPTION_H_

#include <exception>
#include <string>

namespace uu {
namespace core {

/**
 * Exception thrown when there is an attempt to create an object that already exists.
 */
class
    DuplicateElementException: public std::exception
{

  public:

    /**
     * @param value a string indicating type and name (or other identifier) of the element that determined the exception.
     */
    DuplicateElementException(
        std::string value
    );

    ~DuplicateElementException(
    )
    throw ();

    /**
     * Information about the exception.
     * @return an error message describing the occurred problem
     */
    virtual
    const char*
    what(
    ) const
    throw();

  private:

    std::string value;

};

}
}

#endif
