#include "core/exceptions/OperationNotSupportedException.hpp"

namespace uu {
namespace core {

OperationNotSupportedException::OperationNotSupportedException(std::string value)
{
    OperationNotSupportedException::value = "Operation not supported: " + value;
}

OperationNotSupportedException::~OperationNotSupportedException() throw () {}

const char*
OperationNotSupportedException::what() const throw()
{
    return value.data();
}

} // namespace core
} // namespace uu
