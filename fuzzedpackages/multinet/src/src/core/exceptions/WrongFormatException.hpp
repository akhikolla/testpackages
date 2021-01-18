/**
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_EXCEPTIONS_WRONGFORMATEXCEPTION_H_
#define UU_CORE_EXCEPTIONS_WRONGFORMATEXCEPTION_H_

#include <exception>
#include <string>

namespace uu {
namespace core {

/**
 * Exception thrown when a data file is not correctly formatted.
 */
class
    WrongFormatException: public std::exception
{

  public:

    /**
     * @param value message explaining what was wrong in the format
     */
    WrongFormatException(
        std::string value
    );

    ~WrongFormatException(
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
