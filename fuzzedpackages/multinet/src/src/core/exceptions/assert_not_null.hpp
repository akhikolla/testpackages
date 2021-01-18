/**
 * History:
 * - 2018.08.11 file created
 */

#ifndef UU_CORE_EXCEPTIONS_ASSERT_H_
#define UU_CORE_EXCEPTIONS_ASSERT_H_

#include "core/exceptions/NullPtrException.hpp"
#include <string>

namespace uu {
namespace core {

/**
 * Throws a NullPointerException when the first parameter is null.
 * @param ptr the pointer to be checked for NULL
 * @param function the name of the function where the exception is called
 * @param par the name of the function parameter
 */
void
assert_not_null(
    const void* ptr,
    std::string function,
    std::string par
);

}
}

#endif
