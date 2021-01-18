#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace net {

template<typename G>
size_t
order(
    const G* g
)
{
    core::assert_not_null(g, "order", "g");
    return g->vertices()->size();
}


}
}

