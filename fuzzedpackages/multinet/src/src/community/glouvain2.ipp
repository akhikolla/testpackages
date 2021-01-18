#include <unordered_map>
#include <vector>
#include "community/_impl/GMetaNetwork.hpp"
#include "community/louvain.hpp"

namespace uu {
namespace net {

template <typename M>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const typename M::layer_type>>>
glouvain2(
    const M* g,
    double omega
)
{

    auto multilayer_metanetwork = convert(g, omega);

    auto meta = std::move(std::get<0>(multilayer_metanetwork));
    auto mapping = std::get<1>(multilayer_metanetwork);

    std::vector<std::unique_ptr<GMetaNetwork>> passes;

    //auto meta = pass(meta1.get());

    while (meta)
    {
        //std::cout << "pass" << std::endl;
        passes.push_back(std::move(meta));
        auto w = passes.back().get();
        meta = pass(w);
    }


    auto c = communities(passes);

    auto communities = std::make_unique<CommunityStructure<VertexLayerCommunity<const typename M::layer_type>>>();

    for (auto meta_community: *c)
    {
        auto community = std::make_unique<VertexLayerCommunity<const typename M::layer_type>>();

        for (auto meta_vertex: *meta_community)
        {
            community->add(mapping.at(meta_vertex));
        }

        communities->add(std::move(community));
    }

    return communities;
}

}
}
