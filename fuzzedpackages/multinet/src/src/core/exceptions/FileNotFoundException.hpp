/**
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_EXCEPTIONS_FILENOTFOUNDEXCEPTION_H_
#define UU_CORE_EXCEPTIONS_FILENOTFOUNDEXCEPTION_H_

#include <exception>
#include <string>

namespace uu {
namespace core {

/**
 * Exception thrown when a non-existing file is opened.
 */
class
    FileNotFoundException : public std::exception
{

  public:

    /**
     * @param value path of the file
     */
    FileNotFoundException(
        std::string value
    );

    ~FileNotFoundException(
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
