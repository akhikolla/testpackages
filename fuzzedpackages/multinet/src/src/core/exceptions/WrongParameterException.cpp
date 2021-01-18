#include "core/exceptions/WrongParameterException.hpp"

namespace uu {
namespace core {

WrongParameterException::WrongParameterException(std::string value)
{
    WrongParameterException::value = "Wrong parameter: " + value;
}

WrongParameterException::~WrongParameterException() throw() {}

const char*
WrongParameterException::what() const throw()
{
    return value.data();
}

} // namespace core
} // namespace uu
