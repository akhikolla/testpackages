namespace uu {
namespace net {


template <typename M, typename L>
std::unique_ptr<CommunityStructure<PillarCommunity<L>>>
abacus(
    const M* mnet,
    int min_actors,
    int min_layers
)
{
    std::vector<std::unique_ptr<CommunityStructure<Community<const Vertex*>>>> coms;
    std::unordered_map<const L*, CommunityStructure<Community<const Vertex*>>*> single_layer_communities;

    for (auto layer: *mnet->layers())
    {
        auto c = label_propagation(layer);
        single_layer_communities[layer] = c.get();
        coms.push_back(std::move(c));
    }

    return eclat_merge(mnet, single_layer_communities, min_actors, min_layers);
}

}
}

