/**
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_EXCEPTIONS_ELEMENTNOTFOUNDEXCEPTION_H_
#define UU_CORE_EXCEPTIONS_ELEMENTNOTFOUNDEXCEPTION_H_

#include <exception>
#include <string>

namespace uu {
namespace core {

/**
 * Exception thrown when an object required by a function does not exist.
 */
class
    ElementNotFoundException: public std::exception
{

  public:

    /**
     * @param value a string indicating type and name (or other identifier) of the element not found
     */
    ElementNotFoundException(
        std::string value
    );

    ~ElementNotFoundException(
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
