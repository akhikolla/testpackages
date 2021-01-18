#include "networks/OrderedMultiplexNetwork.hpp"

#include "networks/_impl/observers/PropagateObserver.hpp"
#include "networks/_impl/observers/LayerObserver.hpp"
#include "networks/_impl/observers/PropagateAddEraseObserver.hpp"

namespace uu {
namespace net {


OrderedMultiplexNetwork::
OrderedMultiplexNetwork(
    const std::string& name
)
{
    auto vs = std::make_unique<AttrVertexStore>();

    auto ls = std::make_unique<VertexOverlappingOrderedLayerStore<Network>>();

    using EA = Attributes<MLEdge<Vertex,Network>, UserDefinedAttrs<MLEdge<Vertex,Network>>>;
    auto e_attr = std::make_unique<EA>();
    auto es = std::make_unique<AttributedDynamicInterlayerSimpleEdgeStore<Vertex,Network,EA>>(std::move(e_attr));

    // @todo add observers

    TMultilayerNetworkType t;

    data_ = std::make_unique<TMultilayerNetwork<AttrVertexStore, MLOrderedLayerStore, MLSimpleEdgeStore>>(name,t,std::move(vs),std::move(ls),std::move(es));

}


AttrVertexStore*
OrderedMultiplexNetwork::
actors(
)
{
    return data_->vertices();
}


const AttrVertexStore*
OrderedMultiplexNetwork::
actors(
) const
{
    return data_->vertices();
}


MLOrderedLayerStore*
OrderedMultiplexNetwork::
layers(
)
{
    return data_->layers();
}


const MLOrderedLayerStore*
OrderedMultiplexNetwork::
layers(
) const
{
    return data_->layers();
}

std::string
OrderedMultiplexNetwork::
summary(
) const
{

    size_t num_intra_edges = 0;

    for (auto&& layer: *layers())
    {
        num_intra_edges += layer->edges()->size();
    }

    size_t num_inter_edges = 0;

    size_t num_actors = actors()->size();

    size_t num_layers = layers()->size();

    size_t num_nodes = 0;

    for (auto&& layer: *layers())
    {
        num_nodes += layer->vertices()->size();
    }

    size_t num_edges = num_intra_edges + num_inter_edges;

    std::string summary =
        "Multilayer Network [" +
        std::to_string(num_actors) + (num_actors==1?" actor, ":" actors, ") +
        std::to_string(num_layers) + (num_layers==1?" layer, ":" layers, ") +
        std::to_string(num_nodes) + (num_nodes==1?" vertex, ":" vertices, ") +
        std::to_string(num_edges) + (num_edges==1?" edge ":" edges ") +
        "(" + std::to_string(num_intra_edges) + "," +  std::to_string(num_inter_edges) + ")]";
    return summary;
}


}
}

