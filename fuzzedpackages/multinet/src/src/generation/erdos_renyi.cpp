#include "generation/erdos_renyi.hpp"

#include "core/exceptions/OperationNotSupportedException.hpp"
#include "core/utils/random.hpp"
#include "generation/standard_graphs.hpp"
#include "generation/utils.hpp"
#include "operations/_impl/add_predefined_subgraphs.hpp"

namespace uu {
namespace net {

std::unique_ptr<Network>
erdos_renyi_nm(
    size_t n,
    size_t m
)
{
    std::string name = "ER";

    auto g = std::make_unique<Network>(name);
    add_vertices(g.get(), n);

    auto edge_ids = core::get_k_uniform(n*(n-1)/2, m);

    for (auto edge_id:  edge_ids)
    {
        size_t v_id1 = 0;

        while (edge_id >= n - v_id1 - 1)
        {
            edge_id -= n - v_id1 - 1;
            v_id1++;
        }

        size_t v_id2 = edge_id + v_id1 + 1;

        auto v1 = g->vertices()->at(v_id1);
        auto v2 = g->vertices()->at(v_id2);

        g->edges()->add(v1, v2);
    }

    return g;
}



std::unique_ptr<Network>
erdos_renyi_np(
    size_t n,
    double p
)
{
    size_t max_edges = n*(n-1)/2;
    size_t m = core::get_binomial(max_edges, p);
    return erdos_renyi_nm(n, m);

}


}
}

