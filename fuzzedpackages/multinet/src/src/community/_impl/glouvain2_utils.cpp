
#include "community/_impl/glouvain2_utils.hpp"

#include <chrono>
#include <unordered_map>
#include <vector>
#include "community/CommunityStructure.hpp"
#include "community/Community.hpp"
#include "community/_impl/GMetaNetwork.hpp"
#include "objects/EdgeMode.hpp"
#include "objects/Vertex.hpp"
#include "measures/degree.hpp"

namespace uu {
namespace net {


std::tuple<std::unique_ptr<GMetaNetwork>, std::map<const Vertex*, std::pair<const Vertex*, const Network*>>, std::vector<std::unique_ptr<const Vertex>>>
        convert(
            const MultilayerNetwork* g,
            double omega
        )
{
    std::map<std::pair<const Vertex*, const Network*>, const Vertex*> mapping;
    std::map<const Vertex*, std::pair<const Vertex*, const Network*>> reverse_mapping;

    //std::map<std::pair<const Vertex*, const Network*>, size_t> deg;
    //std::map<std::pair<const Vertex*, const Network*>, double> ome;
    //std::map<const Network*, size_t> m;

    auto meta = std::make_unique<GMetaNetwork>();

    std::vector<std::unique_ptr<const Vertex>> metavertices;


    size_t v_id = 0;

    /*
     double mu = 0;

    for (auto l: *g->layers())
    {
        m[l] = size(l);
    }
    for (auto v: *g->vertices())
    {
        for (size_t i = 0; i < g->layers()->size(); i++)
        {
            auto l1 = g->layers()->at(i);
            meta->add();

            if (!l1->vertices()->contains(v))
                continue;
            auto iv1 = std::make_pair(v, l1);
            auto d = degree(l1, v);
            deg[iv1] = d;
            ome[iv1] = 0;
            mu += d/2.0;
        }
    }
    for (auto v: *g->vertices())
    {
        for (auto l1: *g->layers())
        {
            if (!l1->vertices()->contains(v))
                continue;
            auto iv1 = std::make_pair(v, l1);
            for (auto l2: *g->layers())
            {
                if (l1>=l2)
                    continue;
                if (!l2->vertices()->contains(v))
                    continue;
                auto iv2 = std::make_pair(v, l2);

                ome[iv1] += omega;
                ome[iv2] += omega;

                mu += omega;
            }
        }
    }

     */
    for (size_t i = 0; i < g->layers()->size(); i++)
    {
        auto l = g->layers()->at(i);

        for (auto v: *l->vertices())
        {
            auto metavertex = std::make_unique<const Vertex>(std::to_string(v_id++));
            auto intralayer_vertex = std::make_pair(v, l);
            mapping[intralayer_vertex] = metavertex.get();
            reverse_mapping[metavertex.get()] = intralayer_vertex;
            meta->add(metavertex.get());
            metavertices.push_back(std::move(metavertex));
            ////std::cout << "adding " << (*v) << "@" << l->name <<
            //" as " << (*metavertex) << std::endl;
        }

        for (auto e: *l->edges())
        {

            auto iv1 = std::make_pair(e->v1, l);
            auto iv2 = std::make_pair(e->v2, l);
            auto v1 = mapping.at(iv1);
            auto v2 = mapping.at(iv2);

            meta->edge(v1, v2, i+1);

            //std::cout << "adding " << (*v1) << "--" << (*v2) << " on layer " << i << std::endl;
        }
    }

    for (auto v: *g->actors())
    {
        for (auto l1: *g->layers())
        {
            if (!l1->vertices()->contains(v))
            {
                continue;
            }

            for (auto l2: *g->layers())
            {
                if (l1>=l2)
                {
                    continue;
                }

                if (!l2->vertices()->contains(v))
                {
                    continue;
                }

                auto iv1 = std::make_pair(v, l1);
                auto iv2 = std::make_pair(v, l2);
                auto v1 = mapping.at(iv1);
                auto v2 = mapping.at(iv2);

                meta->edge(v1, v2, 0, omega);

                //std::cout << "adding " << (*v1) << "--" << (*v2) << " intralayer" << std::endl;
            }
        }
    }


    ////std::cout << "conversion done!" << std::endl;


    return std::make_tuple(std::move(meta), reverse_mapping, std::move(metavertices));
}


std::tuple<std::unique_ptr<GMetaNetwork>, std::map<const Vertex*, std::pair<const Vertex*, const Network*>>, std::vector<std::unique_ptr<const Vertex>>>
        convert(
            const OrderedMultiplexNetwork* g,
            double omega
        )
{
    std::map<std::pair<const Vertex*, const Network*>, const Vertex*> mapping;
    std::map<const Vertex*, std::pair<const Vertex*, const Network*>> reverse_mapping;

    auto meta = std::make_unique<GMetaNetwork>();

    std::vector<std::unique_ptr<const Vertex>> metavertices;


    size_t v_id = 0;

    for (size_t i = 0; i < g->layers()->size(); i++)
    {
        auto l = g->layers()->at(i);

        for (auto v: *l->vertices())
        {
            auto metavertex = std::make_unique<const Vertex>(std::to_string(v_id++));
            auto intralayer_vertex = std::make_pair(v, l);
            mapping[intralayer_vertex] = metavertex.get();
            reverse_mapping[metavertex.get()] = intralayer_vertex;
            meta->add(metavertex.get());
            metavertices.push_back(std::move(metavertex));
            ////std::cout << "adding " << (*v) << "@" << l->name <<
            //" as " << (*metavertex) << std::endl;
        }

        for (auto e: *l->edges())
        {

            auto iv1 = std::make_pair(e->v1, l);
            auto iv2 = std::make_pair(e->v2, l);
            auto v1 = mapping.at(iv1);
            auto v2 = mapping.at(iv2);

            meta->edge(v1, v2, i+1);

        }
    }

    for (auto v: *g->actors())
    {
        for (size_t i = 0; i < g->layers()->size()-1; i++)
        {
            auto l1 = g->layers()->at(i);
            auto l2 = g->layers()->at(i+1);

            if (!l1->vertices()->contains(v))
            {
                continue;
            }

            if (!l2->vertices()->contains(v))
            {
                continue;
            }

            auto iv1 = std::make_pair(v, l1);
            auto iv2 = std::make_pair(v, l2);
            auto v1 = mapping.at(iv1);
            auto v2 = mapping.at(iv2);

            meta->edge(v1, v2, 0, omega);
        }
    }

    return std::make_tuple(std::move(meta), reverse_mapping, std::move(metavertices));
}


void
expand(
    const std::vector<std::unique_ptr<GMetaNetwork>>& levels,
    size_t i,
    const Vertex* v,
    Community<const Vertex*>* com
)
{
    if (i==0)
    {
        for (auto original_vertex: levels.at(i)->mapping.at(v))
        {
            ////std::cout << "adding " << (*original_vertex) << std::endl;
            com->add(original_vertex);
        }
    }

    else
    {
        for (auto previous_vertex: levels.at(i)->mapping.at(v))
        {
            expand(levels, i-1, previous_vertex, com);
        }
    }
}

std::unique_ptr<CommunityStructure<Community<const Vertex*>>>
communities(
    const std::vector<std::unique_ptr<GMetaNetwork>>& levels
)
{
    auto res = std::make_unique<CommunityStructure<Community<const Vertex*>>>();

    size_t i = levels.size()-1;

    for (auto v: *levels.at(i)->get()->vertices())
    {
        auto community = std::make_unique<Community<const Vertex*>>();

        expand(levels, i, v, community.get());

        res->add(std::move(community));
    }

    return res;
}


std::unique_ptr<GMetaNetwork>
aggregate(
    const GMetaNetwork* meta,
    std::unordered_map<const Vertex*, size_t> community
)
{
    auto meta_agg = std::make_unique<GMetaNetwork>();

    std::unordered_map<size_t, std::set<const Vertex*>> vertices;

    for (auto pair: community)
    {
        vertices[pair.second].insert(pair.first);
        /*
         auto it = vertices.find(pair->second);
         if (it == vertices.end())
         {
         improvement[neighbor_community] += contribution;
         }
         else {
         it->second += contribution;
         }
         */
    }


    for (auto pair: vertices)
    {
        meta_agg->add(pair.second.begin(), pair.second.end());

        /*for (auto v: pair.second)
         {
         //std::cout << (*v) << " -> " << (*meta_v) << std::endl;
         }*/
    }

    for (auto e: *meta->get()->edges())
    {
        auto type = meta->get_type(e);
        double weight = meta->get_weight(e);
        meta_agg->edge(e->v1, e->v2, type, weight);

    }

    return meta_agg;
}


std::unique_ptr<GMetaNetwork>
pass(
    const GMetaNetwork* meta
)
{
    //std::cout << "PASS" << std::endl;

    std::unordered_map<const Vertex*, size_t> community;
    std::unordered_map<size_t, std::set<const Vertex*>> vertices_in_community;
    std::map<std::pair<const Vertex*, size_t>, double> w_degree;
    size_t comm_id = 0;

    // initializing communities
    for (auto v: *meta->get()->vertices())
    {
        ////std::cout << (*v) << ": " << comm_id << std::endl;
        vertices_in_community[comm_id].insert(v);
        community[v] = comm_id;
        comm_id++;
    }

    //double mu = 0;
    // compute the weight of each layer and strength of each vertex
    std::unordered_map<int, double> m;

    for (auto e: *meta->get()->edges())
    {
        auto edge_type = meta->get_type(e);
        auto edge_weight = meta->get_weight(e);

        auto vertex_layer1 = std::make_pair(e->v1, edge_type);
        w_degree[vertex_layer1] += edge_weight;

        auto vertex_layer2 = std::make_pair(e->v2, edge_type);
        w_degree[vertex_layer2] += edge_weight;

        m[edge_type] += edge_weight;

        ////std::cout << " m " << m << " " << g->get_weight(e).null << std::endl;
    }

    bool change = false;
    bool improved = false;

    do
    {
        change = false;

        for (auto v: *meta->get()->vertices())
        {

            //std::cout << "Vertex " << (*v) << ":" << std::endl;

            auto current_community = community.at(v);
            // communities of vertices adjacent to v
            std::set<size_t> neighboring_communities;
            // layers where v has positive degree (0 = interlayer)
            std::set<size_t> positive_degree;

            for (auto e: *meta->get()->edges()->incident(v))
            {
                auto edge_type = meta->get_type(e);
                positive_degree.insert(edge_type);

                auto n = e->v1;

                if (v == n)
                {
                    n = e->v2;
                }

                auto c = community.at(n);

                if (c != current_community)
                {
                    neighboring_communities.insert(c);
                }
            }

            if (neighboring_communities.size() == 0)
            {
                continue;
            }

            std::unordered_map<size_t, double> improvement;

            for (auto c: neighboring_communities)
            {
                improvement[c] = 0;

                for (auto n: vertices_in_community.at(c))
                {

                    //std::cout << (*n) << " (";

                    auto edges = meta->get()->edges()->get(v, n);
                    std::unordered_map<size_t, const Edge*> edge_map;

                    for (auto e: edges)
                    {
                        auto edge_type = meta->get_type(e);
                        edge_map[edge_type] = e;
                    }

                    for (auto edge_type: positive_degree)
                    {
                        double contribution;

                        auto pair = edge_map.find(edge_type);

                        if (edge_type == 0) // interlayer link
                        {
                            if (pair == edge_map.end())
                            {
                                // no interlayer edges between these vertices
                                contribution = 0;
                            }

                            else
                            {
                                auto e = pair->second;
                                contribution = meta->get_weight(e);
                            }
                        }

                        else
                        {
                            double A_ij;

                            if (pair == edge_map.end())
                            {
                                // no edges between these vertices
                                A_ij = 0;
                            }

                            else
                            {
                                auto e = pair->second;
                                A_ij = meta->get_weight(e);
                            }

                            auto vertex_layer1 = std::make_pair(v, edge_type);
                            auto vertex_layer2 = std::make_pair(n, edge_type);

                            auto d1 = w_degree[vertex_layer1];
                            auto d2 = w_degree[vertex_layer2];

                            contribution = A_ij - d1*d2/m[edge_type]/2;

                        }

                        //std::cout << contribution << ") ";
                        improvement[c] += contribution;
                    }
                }

                //std::cout << ": " << improvement[c] << std::endl;
            }


            double loss = 0; // price for leaving current community

            for (auto n: vertices_in_community.at(current_community))
            {
                if (v == n)
                {
                    continue;
                }

                //std::cout << (*n) << " (";

                auto edges = meta->get()->edges()->get(v, n);
                std::unordered_map<size_t, const Edge*> edge_map;

                for (auto e: edges)
                {
                    auto edge_type = meta->get_type(e);
                    edge_map[edge_type] = e;
                }

                for (auto edge_type: positive_degree)
                {
                    double contribution;

                    auto pair = edge_map.find(edge_type);

                    if (edge_type == 0) // interlayer link
                    {
                        if (pair == edge_map.end())
                        {
                            // no interlayer edges between these vertices
                            contribution = 0;
                        }

                        else
                        {
                            auto e = pair->second;
                            contribution = meta->get_weight(e);
                        }
                    }

                    else
                    {
                        double A_ij;

                        if (pair == edge_map.end())
                        {
                            // no edges between these vertices
                            A_ij = 0;
                        }

                        else
                        {
                            auto e = pair->second;
                            A_ij = meta->get_weight(e);
                        }

                        auto vertex_layer1 = std::make_pair(v, edge_type);
                        auto vertex_layer2 = std::make_pair(n, edge_type);

                        auto d1 = w_degree[vertex_layer1];
                        auto d2 = w_degree[vertex_layer2];

                        contribution = A_ij - d1*d2/m[edge_type]/2;

                    }

                    //std::cout << contribution << ") ";
                    loss += contribution;
                }
            }


            //std::cout << "loss: " << loss << std::endl;


            // find best choice
            auto new_community = current_community;
            double current_improvement = 0;

            for (auto pair: improvement)
            {
                if (pair.second - loss > current_improvement)
                {
                    new_community = pair.first;
                    current_improvement = pair.second - loss;
                }
            }

            if (new_community != current_community)
            {
                change = true;
                improved = true;
                vertices_in_community.at(current_community).erase(v);
                vertices_in_community.at(new_community).insert(v);
                community.at(v) = new_community;

                ////std::cout << (*v) << " from " << current_community << " to " << new_community << std::endl;
            }
        }

        //std::cout << "---" << std::endl;
    }
    while (change);

    /*for (auto pair: community)
     {
     //std::cout << (*pair.first) << ": " << pair.second << std::endl;
     }*/

    if (!improved)
    {
        return nullptr;
    }

    auto meta_agg = aggregate(meta, community);

    return meta_agg;
}

}
}

