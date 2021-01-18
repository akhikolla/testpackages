#include "core/exceptions/OutOfBoundsException.hpp"

namespace uu {
namespace core {

OutOfBoundsException::OutOfBoundsException(std::string value)
{
    OutOfBoundsException::value = "Requested element out of bounds: " + value;
}

OutOfBoundsException::OutOfBoundsException() throw () {}

const char*
OutOfBoundsException::what() const throw()
{
    return value.data();
}

} // namespace core
} // namespace uu
