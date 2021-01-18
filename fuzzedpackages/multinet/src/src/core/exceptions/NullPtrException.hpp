/**
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_EXCEPTIONS_NULLPTREXCEPTION_H_
#define UU_CORE_EXCEPTIONS_NULLPTREXCEPTION_H_

#include <exception>
#include <string>

namespace uu {
namespace core {

/**
 * Exception thrown when a pointer parameter passed to a function is null.
 */
class
    NullPtrException: public std::exception
{

  public:

    /**
     * @param value a string indicating the type of element whose pointer is null
     */
    NullPtrException(
        std::string value
    );

    ~NullPtrException(
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
