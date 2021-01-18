/**
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_EXCEPTIONS_WRONGPARAMETEREXCEPTION_H_
#define UU_CORE_EXCEPTIONS_WRONGPARAMETEREXCEPTION_H_

#include <exception>
#include <string>

namespace uu {
namespace core {

/**
 * Exception thrown when a function is called with a parameter value that is not
 * among the ones recognized by the function.
 */
class
    WrongParameterException: public std::exception
{

  public:

    /**
     * @param value a string indicating the name and value of the parameter
     */
    WrongParameterException(
        std::string value
    );

    ~WrongParameterException(
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
