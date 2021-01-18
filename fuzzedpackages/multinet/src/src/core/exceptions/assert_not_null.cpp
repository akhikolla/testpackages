#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace core {

void
assert_not_null(
    const void* ptr,
    std::string function,
    std::string param
)
{
    if (!ptr)
    {
        std::string msg = "function " + function +
                          ", parameter " + param;
        throw NullPtrException(msg);
    }
}

}
}
