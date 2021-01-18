#include <vector>
#include <algorithm>
#include "core/exceptions/assert_not_null.hpp"
#include "core/exceptions/WrongParameterException.hpp"

namespace uu {
namespace net {

template<typename G>
double
exp_degree(
    const G* g,
    const Vertex* v,
    const EdgeMode mode
)
{
    core::assert_not_null(g, "exp_degree", "g");
    core::assert_not_null(g, "exp_degree", "v");

    if (!g->is_probabilistic())
    {
        throw core::WrongParameterException("expected degree can only be computed on probabilistic networks");
    }

    double exp_d = 0;

    for (auto edge: *g->edges()->incident(v, mode))
    {
        auto p = g->get_prob(edge);

        if (!p.null)
        {
            exp_d += p.value;
        }
    }

    return exp_d;
}


}
}

