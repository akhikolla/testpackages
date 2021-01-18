#ifndef UU_COMMUNITY_MLCPM_H_
#define UU_COMMUNITY_MLCPM_H_

#include <memory>
#include <unordered_set>
#include <vector>
#include "core/utils/math.hpp"
#include "core/utils/pretty_printing.hpp"
#include "community/MLCPMCommunity.hpp"
#include "community/CommunityStructure.hpp"

namespace uu {
namespace net {

template <typename M>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const typename M::layer_type>>>
mlcpm(
    const M* mnet,
    size_t k,
    size_t m
);

/**
 * This method finds all the maximal cliques (1) being a superset of clique A
 * (2) extended with actors in B (3) not containing actors in C (4) on
 * at least k actors and m layers.
 */
template <typename M>
std::unordered_set<std::shared_ptr<MultiplexClique<M>>>
find_max_cliques(
    const M* mnet,
    size_t k,
    size_t m
);


template <typename M>
std::unordered_set<const typename M::layer_type*>
neighboring_layers(
    const M* mnet,
    const Vertex* actor1,
    const Vertex* actor2
);

template <typename M>
std::map<std::shared_ptr<MultiplexClique<M>>,std::unordered_set<std::shared_ptr<MultiplexClique<M>>> >
        build_max_adjacency_graph(
            const std::unordered_set<std::shared_ptr<MultiplexClique<M>>>& C,
            size_t k,
            size_t m
        );

template <typename M>
std::unordered_set<std::shared_ptr<MLCPMCommunity<M>>>
find_max_communities(
    //const M* mnet,
    const std::map<std::shared_ptr<MultiplexClique<M>>,std::unordered_set<std::shared_ptr<MultiplexClique<M>>> >& adjacency,
    size_t m
);


template <typename M>
void
find_max_communities(
    const std::map<std::shared_ptr<MultiplexClique<M>>,std::unordered_set<std::shared_ptr<MultiplexClique<M>>> >& adjacency,
    std::shared_ptr<MLCPMCommunity<M>>& A,
    std::vector<std::shared_ptr<MultiplexClique<M>>> Candidates,
    std::unordered_set<std::shared_ptr<MultiplexClique<M>>>& processedCliques,
    layer_sets<M>& processedLayerCombinations,
    size_t m,
    std::unordered_set<std::shared_ptr<MLCPMCommunity<M>>>& result
);


}
}

#include "ml-cpm.ipp"

#endif
