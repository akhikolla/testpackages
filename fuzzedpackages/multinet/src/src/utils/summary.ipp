#include "networks/MultilayerNetwork.hpp"
#include "core/exceptions/OperationNotSupportedException.hpp"
#include "core/exceptions/assert_not_null.hpp"
#include "measures/order.hpp"
#include "measures/size.hpp"

namespace uu {
namespace net {

template<typename G>
std::string
summary_short(
    const G* g
)
{
    core::assert_not_null(g, "summary_short", "g");

    auto n = std::to_string(order(g));
    auto m = std::to_string(size(g));
    auto dir = g->is_directed()?"Dir":"Und";
    auto loops = g->allows_loops()?", Loops":"";
    return "net(" + n + ", " + m  + ", " + dir + loops + ")";
}

template<>
std::string
summary_short(
    const MultilayerNetwork* net
)
{
    core::assert_not_null(net, "summary_short", "g");

    size_t num_intra_edges = 0;

    for (auto layer: *net->layers())
    {
        num_intra_edges += layer->edges()->size();
    }

    size_t num_inter_edges = net->interlayer_edges()->size();

    size_t num_actors = net->vertices()->size();

    size_t num_layers = net->layers()->size();

    size_t num_vertices = 0;

    for (auto layer: *net->layers())
    {
        num_vertices += layer->vertices()->size();
    }

    size_t num_edges = num_intra_edges + num_inter_edges;

    std::string summary =
        "ml-net[" +
        std::to_string(num_actors) + ", " +
        std::to_string(num_layers) + ", " +
        std::to_string(num_vertices) + ", " +
        std::to_string(num_edges) +
        " (" + std::to_string(num_intra_edges) + "," +  std::to_string(num_inter_edges) + ")]";

    return summary;
}

}
}

