#include "core/attributes/Time.hpp"
#include "core/attributes/conversion.hpp"

namespace uu {
namespace core {


std::ostream&
operator<<(
    std::ostream& os,
    const Time& t
)
{
    os << to_string(t, kDEFAULT_TIME_FORMAT);
    return os;
}


}
}

