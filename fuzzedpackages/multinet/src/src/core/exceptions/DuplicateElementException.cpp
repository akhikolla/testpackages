#include "core/exceptions/DuplicateElementException.hpp"

namespace uu {
namespace core {

DuplicateElementException::DuplicateElementException(std::string value)
{
    DuplicateElementException::value = "Duplicate element: " + value;
}

DuplicateElementException::~DuplicateElementException() throw () {}

const char*
DuplicateElementException::what() const throw()
{
    return value.data();
}

} // namespace core
} // namespace uu
