#include "networks/MultilayerNetwork.hpp"

#include "networks/_impl/observers/PropagateObserver.hpp"
#include "networks/_impl/observers/LayerObserver.hpp"
//#include "networks/_impl/observers/PropagateAddEraseObserver.hpp"

namespace uu {
namespace net {


MultilayerNetwork::
MultilayerNetwork(
    const std::string& name
)
{
    auto vs = std::make_unique<AttrVertexStore>();

    auto ls = std::make_unique<VertexOverlappingLayerStore<Network>>();

    using EA = Attributes<MLEdge<Vertex,Network>, UserDefinedAttrs<MLEdge<Vertex,Network>>>;
    auto e_attr = std::make_unique<EA>();
    auto es = std::make_unique<AttributedDynamicInterlayerSimpleEdgeStore<Vertex,Network,EA>>(std::move(e_attr));

    // @todo missing observer?

    // register an observer to propagate the removal of vertices to the layers
    auto obs1 = std::make_unique<PropagateObserver<VertexOverlappingLayerStore<Network>, const Vertex>>(ls.get());
    vs->attach(obs1.get());

    // register an observer to react to the addition/removal of layers
    auto obs2 = std::make_unique<LayerObserver<AttributedDynamicInterlayerSimpleEdgeStore<Vertex,Network,EA>, Network>>(es.get());
    ls->attach(obs2.get());

    TMultilayerNetworkType t;

    data_ = std::make_unique<TMultilayerNetwork<AttrVertexStore, MLLayerStore, MLSimpleEdgeStore>>(name, t, std::move(vs), std::move(ls), std::move(es));

    data_->register_observer(std::move(obs1));
    data_->register_observer(std::move(obs2));

}


AttrVertexStore*
MultilayerNetwork::
actors(
)
{
    return data_->vertices();
}


const AttrVertexStore*
MultilayerNetwork::
actors(
) const
{
    return data_->vertices();
}


MLLayerStore*
MultilayerNetwork::
layers(
)
{
    return data_->layers();
}


const MLLayerStore*
MultilayerNetwork::
layers(
) const
{
    return data_->layers();
}


MLSimpleEdgeStore*
MultilayerNetwork::
interlayer_edges(
)
{
    return data_->interlayer_edges();
}


const MLSimpleEdgeStore*
MultilayerNetwork::
interlayer_edges(
) const
{
    return data_->interlayer_edges();
}

std::string
MultilayerNetwork::
summary(
) const
{

    size_t num_intra_edges = 0;

    for (auto layer: *layers())
    {
        num_intra_edges += layer->edges()->size();
    }

    size_t num_inter_edges = interlayer_edges()->size();

    size_t num_actors = actors()->size();

    size_t num_layers = layers()->size();

    size_t num_nodes = 0;

    for (auto layer: *layers())
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


bool
MultilayerNetwork::
is_ordered(
) const
{
    return false;
}

bool
MultilayerNetwork::
allows_interlayer_edges(
) const
{
    return false;
}

}
}

