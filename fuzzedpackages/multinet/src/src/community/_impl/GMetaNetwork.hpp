#ifndef UU_NET_COMMUNITY_GMETANETWORK_H_
#define UU_NET_COMMUNITY_GMETANETWORK_H_

#include <unordered_map>
#include <vector>
#include <memory>
#include "networks/MultiNetwork.hpp"
#include "objects/EdgeMode.hpp"
#include "objects/Vertex.hpp"
#include "objects/Edge.hpp"

namespace uu {
namespace net {

class GMetaNetwork
{
  public:

    std::unique_ptr<MultiNetwork> w;
    std::unordered_map<const Edge*, size_t> edge_type;
    std::unordered_map<const Edge*, double> edge_weight;
    size_t order = 0;
    std::unordered_map<const Vertex*, std::set<const Vertex*>> mapping;
    std::unordered_map<const Vertex*, const Vertex*> reverse_mapping;

    GMetaNetwork();

    template <typename VertexIterator>
    const Vertex*
    add(
        VertexIterator begin,
        VertexIterator end
    );


    const Vertex*
    add(
        const Vertex* u
    );


    /* type: -1 (interlayer), n (layer n) */
    const Edge*
    edge(
        const Vertex* u,
        const Vertex* v,
        size_t type,
        double weight = 1.0
    );

    double
    get_weight(
        const Edge*
    ) const;

    size_t
    get_type(
        const Edge*
    ) const;

    const MultiNetwork*
    get(
    ) const;

};

}
}

#include "GMetaNetwork.ipp"

#endif
