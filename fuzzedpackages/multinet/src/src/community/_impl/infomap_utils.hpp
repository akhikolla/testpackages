#ifndef UU_MNET_COMMUNITY_INFOMAPUTILS_H_
#define UU_MNET_COMMUNITY_INFOMAPUTILS_H_

#ifndef NS_INFOMAP
#define NS_INFOMAP
#endif

#include <cstddef> // to prevent '::max_align_t' has not been declared error
#include <cstdio>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include "io/ProgramInterface.h"
#include "io/convert.h"
#include "utils/FileURI.h"
#include "infomap/MultiplexNetwork.h"
#include "io/HierarchicalNetwork.h"
#include "io/Config.h"
#include "core/exceptions/ExternalLibException.hpp"
#include "core/attributes/conversion.hpp"
#include "core/utils/CSVReader.hpp"
#include "objects/Vertex.hpp"
#include "community/VertexLayerCommunity.hpp"
#include "community/CommunityStructure.hpp"

namespace uu {
namespace net {



// Writes the network in the format used by the infomap library
// into filename_prefix+"_uu_net"
template <typename M>
void
multinet_to_infomap(
    const M* net,
    const std::string& filename_prefix
);


template <typename M>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const typename M::layer_type>>>
read_infomap_communities(
    const M* net,
    const std::string& filename_prefix
);



template <typename M>
void
multinet_to_infomap(
    const M* net,
    infomap::MultiplexNetwork& network
)
{

    std::unordered_map<const Vertex*, size_t> actors_ids;
    size_t a_id = 0;

    for (auto a: *net->actors())
    {
        actors_ids[a] = a_id++;
    }

    size_t l_id = 0;

    for (auto l: *net->layers())
    {
        for (auto e: *l->edges())
        {
            network.addMultiplexLink(l_id, actors_ids[e->v1], l_id, actors_ids[e->v2], 1.0);
        }

        l_id++;
    }

    network.finalizeAndCheckNetwork(true);
}

template <typename M>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const typename M::layer_type>>>
to_communities(
    const M* net,
    infomap::HierarchicalNetwork& resultNetwork
)
{
    std::unordered_map<size_t, const Vertex*> actors;
    size_t a_id = 0;

    for (auto a: *net->actors())
    {
        actors[a_id++] = a;
    }

    std::unordered_map<size_t, std::vector<std::pair<const Vertex*,const typename M::layer_type*>>> comm;

    /*
    for (infomap::LeafIterator leafIt(&resultNetwork.getRootNode()); !leafIt.isEnd(); ++leafIt)
    {
        auto actor = actors[leafIt->originalLeafIndex];
        auto c_id = leafIt->parentNode->parentIndex;
        for (auto l: *net->layers())
        {
            if (l->vertices()->contains(actor))
            {
                comm[c_id].push_back(std::make_pair(actor,l));
            }
        }
    }*/
    for (infomap::TreeIterator it(&resultNetwork.getRootNode(), 1); !it.isEnd(); ++it)
    {
        infomap::SNode &node = *it;

        if (node.isLeafNode())
        {
            auto actor = actors[node.originalLeafIndex];
            auto c_id = it.moduleIndex() + 1;

            for (auto l: *net->layers())
            {
                if (l->vertices()->contains(actor))
                {
                    comm[c_id].push_back(std::make_pair(actor,l));
                }
            }
        }
    }


    auto communities = std::make_unique<CommunityStructure<VertexLayerCommunity<const typename M::layer_type>>>();

    for (auto pair: comm)
    {
        auto c = std::make_unique<VertexLayerCommunity<const typename M::layer_type>>();

        for (auto vertex_layer_pair: pair.second)
        {
            c->add(vertex_layer_pair);
        }

        communities->add(std::move(c));
    }

    return communities;
}


template <typename M>
void
multinet_to_infomap(
    const M* net,
    const std::string& filename
)
{

    size_t num_actors = net->actors()->size();
    std::unordered_map<const Vertex*, size_t> actors_ids(num_actors);
    std::vector<const Vertex*> actors;
    size_t a_id = 1;

    for (auto a: *net->actors())
    {
        actors_ids[a] = a_id++;
        actors.push_back(a);
    }

    std::ofstream outfile;
    outfile.open(filename);

    outfile << "# " << filename << std::endl;
    outfile << "*Vertices " << net->actors()->size() << std::endl;

    for (size_t a_id=0; a_id<num_actors; a_id++)
    {
        outfile << (a_id+1) << " \"" << actors.at(a_id)->name << "\"" << std::endl;
    }

    outfile << "*Multilayer" << std::endl;
    outfile << "# layer node layer node [weight]" << std::endl;

    size_t l_id = 1;

    for (auto l: *net->layers())
    {
        for (auto e: *l->edges())
        {
            outfile << l_id << " " << actors_ids[e->v1] << " " << l_id << " " << actors_ids[e->v2] << std::endl;
        }

        l_id++;
    }

    outfile.close();


    /*
    outfile.open(filename_prefix+"_uu_map");
    for (auto p: actors) {
        outfile << p.second.first << " " << p.second.second <<  std::endl;
    }
    outfile.close();
     */
}


template <typename M>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const typename M::layer_type>>>
read_infomap_communities(
    const M* net,
    const std::string& filename
)
{

    // read id -> name mapping for actors
    /*
    std::map<std::string,std::string> actors;

    CSVReader ids;
    ids.setFieldSeparator(' ');
    ids.open(filename_prefix+"_uu_map");
    while (ids.has_next()) {
        std::vector<std::string> f = ids.get_next();
        actors[f[0]] = f[1];
    }
    ids.close();
     */

    std::unordered_map<size_t, const Vertex*> actors;
    size_t a_id = 1;

    for (auto a: *net->actors())
    {
        actors[a_id++] = a;
    }

    std::unordered_map<size_t, std::vector<std::pair<const Vertex*,const typename M::layer_type*>>> comm;

    core::CSVReader clusters;
    clusters.set_field_separator(' ');
    clusters.open(filename);
    clusters.get_next();
    clusters.get_next();

    while (clusters.has_next())
    {
        std::vector<std::string> f = clusters.get_next();
        auto actor = actors[core::to_int(f[0])];
        auto c_id = core::to_int(f[1]);

        for (auto l: *net->layers())
        {
            if (l->vertices()->contains(actor))
            {
                comm[c_id].push_back(std::make_pair(actor,l));
            }
        }
    }

    clusters.close();


    auto communities = std::make_unique<CommunityStructure<VertexLayerCommunity<const typename M::layer_type>>>();

    for (auto pair: comm)
    {
        auto c = std::make_unique<VertexLayerCommunity<const typename M::layer_type>>();

        for (auto vertex_layer_pair: pair.second)
        {
            c->add(vertex_layer_pair);
        }

        communities->add(std::move(c));
    }

    return communities;
}

}
}

#endif
