#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace net {

template<typename G>
size_t
size(
    const G* g
)
{
    core::assert_not_null(g, "size", "g");
    return g->edges()->size();
}


}
}

