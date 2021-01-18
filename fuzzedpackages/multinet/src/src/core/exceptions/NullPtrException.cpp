#include "core/exceptions/NullPtrException.hpp"

namespace uu {
namespace core {

NullPtrException::NullPtrException(std::string value)
{
    NullPtrException::value = "NULL pointer: " + value;
}

NullPtrException::~NullPtrException() throw () {}

const char*
NullPtrException::what() const throw()
{
    return value.data();
}

}
}
