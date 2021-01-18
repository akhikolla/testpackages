#ifndef UU_MNET_COMMUNITY_MODULARITY_H_
#define UU_MNET_COMMUNITY_MODULARITY_H_

#include <unordered_map>
#include <numeric>
#include <utility>
#include "objects/EdgeMode.hpp"

#undef CS

namespace uu {
namespace net {


template <typename M, typename CS>
double
modularity(const M* mnet, const CS* communities, double omega)
{

    double res = 0;
    double mu = 0;

    std::unordered_map<std::string,int> m_s;

    for (auto s: *mnet->layers())
    {
        double m = s->edges()->size();

        if (!s->is_directed())
        {
            m *= 2;
        }

        // FIX TO THE ORIGINAL EQUATION WHEN THERE ARE NO EDGES
        if (m == 0)
        {
            m = 1;    // no effect on the formula
        }

        m_s[s->name] = m;
        mu += m;
    }

    for (auto community: *communities)
    {
        for (auto i: *community)
        {
            for (auto j: *community)
            {
                if (i==j)
                {
                    continue;    // not in the original definition - we do this assuming to deal with simple graphs
                }

                //std::cout << i->to_string() << " " << groups.count(i) << std::endl;
                //std::cout << j->to_string() << " " << groups.count(j) << std::endl;

                // same layer
                if (i.second==j.second)
                {
                    //std::cout << "Same group!" << std::endl;
                    //if (mnet.getNetwork(net)->containsEdge(*v_i,*v_j))
                    //    std::cout << "Edge" << std::endl;
                    long k_i = i.second->edges()->neighbors(i.first,EdgeMode::OUT)->size();
                    long k_j = j.second->edges()->neighbors(j.first,EdgeMode::IN)->size();
                    int a_ij = i.second->edges()->get(i.first, j.first)? 1.0 : 0.0;
                    res += a_ij - (double)k_i * k_j / (m_s.at(i.second->name));
                    //std::cout << i->actor->name << " " << j->actor->name << " " << i->layer->name << " "<< k_i << " " <<  k_j << " " <<  m_s.at(i->layer) << std::endl;
                    //std::cout << "->" << res << std::endl;
                }

                if (i.first==j.first)
                {
                    res += omega;
                }
            }
        }

        //std::cout << "->" << m_net << std::endl;
    }

    //std::cout << "same" << std::endl;

    //std::cout << mu << std::endl;
    int num_layers = mnet->layers()->size();

    for (auto actor: *mnet->actors())
    {
        mu+=num_layers*(num_layers-1)*omega; // unclear if we should multiply by c
    }

    //std::cout << mu << std::endl;

    return 1 / mu * res;
}



template <typename M, typename CS>
double
ordered_modularity(const M* mnet, const CS* communities, double omega)
{



    double res = 0;
    double mu = 0;

    std::unordered_map<std::string,int> m_s;

    for (auto s: *mnet->layers())
    {
        double m = s->edges()->size();

        if (!s->is_directed())
        {
            m *= 2;
        }

        // FIX TO THE ORIGINAL EQUATION WHEN THERE ARE NO EDGES
        if (m == 0)
        {
            m = 1;    // no effect on the formula
        }

        m_s[s->name] = m;
        mu += m;
    }

    // quick fix
    std::map<const uu::net::Network*,size_t> layer_index; // FIX ME!!

    for (size_t l_idx=0; l_idx<mnet->layers()->size(); l_idx++)
    {
        layer_index[mnet->layers()->at(l_idx)] = l_idx;
    }

    for (auto community: *communities)
    {
        for (auto i: *community)
        {
            for (auto j: *community)
            {
                if (i==j)
                {
                    continue;    // not in the original definition - we do this assuming to deal with simple graphs
                }

                //std::cout << i->to_string() << " " << groups.count(i) << std::endl;
                //std::cout << j->to_string() << " " << groups.count(j) << std::endl;

                // same layer
                if (i.second==j.second)
                {
                    //std::cout << "Same group!" << std::endl;
                    //if (mnet.getNetwork(net)->containsEdge(*v_i,*v_j))
                    //    std::cout << "Edge" << std::endl;
                    long k_i = i.second->edges()->neighbors(i.first,EdgeMode::OUT)->size();
                    long k_j = j.second->edges()->neighbors(j.first,EdgeMode::IN)->size();
                    int a_ij = i.second->edges()->get(i.first, j.first)? 1.0 : 0.0;
                    res += a_ij - (double)k_i * k_j / (m_s.at(i.second->name));
                    //std::cout << i->actor->name << " " << j->actor->name << " " << i->layer->name << " "<< k_i << " " <<  k_j << " " <<  m_s.at(i->layer) << std::endl;
                    //std::cout << "->" << res << std::endl;
                }

                if (i.first==j.first)
                {
                    if (std::abs((int)(layer_index[i.second]-layer_index[j.second]))<2)
                    {
                        res += omega;
                    }
                }
            }
        }

        //std::cout << "->" << m_net << std::endl;
    }

    //std::cout << "same" << std::endl;

    //std::cout << mu << std::endl;
    int num_layers = mnet->layers()->size();
    int num_actors = mnet->actors()->size();

    mu += 2*(num_layers-1)*num_actors*omega; // unclear if we should multiply by c


    //std::cout << mu << std::endl;

    return 1 / mu * res;
}


}
}


#endif
