#include "operations/slice.hpp"
#include "networks/Network.hpp"
#include "objects/EdgeDir.hpp"
#include "core/attributes/conversion.hpp"
#include "core/attributes/Value.hpp"
#include "core/attributes/Time.hpp"
#include "core/exceptions/assert_not_null.hpp"
#include "core/exceptions/OperationNotSupportedException.hpp"
#include <vector>
#include <fstream>
#include <algorithm>

namespace uu {
namespace net {


std::unique_ptr<OrderedMultiplexNetwork>
slice_equal_time(
    const TemporalNetwork* tnet,
    size_t num_partitions
)
{

    core::assert_not_null(tnet, "slice_equal_time", "tnet");

    auto mpx = std::make_unique<OrderedMultiplexNetwork>(tnet->name);
    std::vector<const Edge*> sorted_edge_vector;
    std::vector<std::vector<const Edge*>> partitioned_edge_vector;

    // create ordered layers for each time partition
    for (size_t i = 0; i<num_partitions; i++)
    {
        EdgeDir dir = (tnet->is_directed()? EdgeDir::DIRECTED : EdgeDir::UNDIRECTED);
        auto g = std::make_unique<uu::net::Network>("l" + std::to_string(i), dir);
        mpx->layers()->push_back(std::move(g));
    }

    // adding all vertices to all the layers
    for (auto vertex: *tnet->vertices())
    {
        mpx->actors()->add(vertex);

        for (auto layer: *mpx->layers())
        {
            layer->vertices()->add(vertex);
        }
    }

    auto max_time = tnet->get_max_time().value;
    auto min_time = tnet->get_min_time().value;

    auto split_time = (max_time - min_time) / (float)num_partitions;

    if (max_time == min_time)
    {
        throw core::OperationNotSupportedException("cannot slice a network with no temporal extension");
    }

    for (auto e : *tnet->edges())
    {
        auto t = tnet->get_time(e);

        if (t.null)
        {
            continue;
        }

        size_t idx = ((t.value - min_time) / split_time);

        if (idx == num_partitions)
        {
            idx--;
        }

        auto layer = mpx->layers()->at(idx);

        layer->edges()->add(e->v1, e->v2);
    }

    return mpx;
}

}
}

