#ifndef UU_MEASURES_LAYER_H_
#define UU_MEASURES_LAYER_H_

#include <vector>
#include <algorithm>
#include <unordered_set>
#include "core/exceptions/assert_not_null.hpp"
#include "core/datastructures/propertymatrix/PropertyMatrix.hpp"
#include "objects/Vertex.hpp"
#include "objects/Dyad.hpp"
#include "objects/Triad.hpp"
#include "objects/EdgeMode.hpp"

namespace uu {
namespace net {

template <typename M>
core::PropertyMatrix<const Vertex*,const typename M::layer_type*,bool>
actor_existence_property_matrix(
    const M* mnet
);

template <typename M>
core::PropertyMatrix<std::pair<const typename M::vertex_type*,const typename M::vertex_type*>,const typename M::layer_type*,bool>
edge_existence_property_matrix(
    const M* mnet
);

template <typename M>
core::PropertyMatrix<Triad,const typename M::layer_type*,bool>
triangle_existence_property_matrix(
    const M* mnet
);

template <typename M>
core::PropertyMatrix<const Vertex*,const typename M::layer_type*,double>
actor_degree_property_matrix(
    const M* mnet, EdgeMode mode
);

template <typename M>
core::PropertyMatrix<const Vertex*,const typename M::layer_type*,double>
actor_cc_property_matrix(
    const M* mnet
);



template <typename M>
core::PropertyMatrix<const Vertex*,const typename M::layer_type*,bool>
actor_existence_property_matrix(
    const M* mnet
)
{
    core::PropertyMatrix<const Vertex*,const typename M::layer_type*,bool> P(mnet->actors()->size(),mnet->layers()->size(),false);

    for (auto layer: *mnet->layers())
        for (auto actor: *layer->vertices())
        {
            P.set(actor,layer,true);
        }

    return P;
}


// @todo document: does not consider self edges
template <typename M>
core::PropertyMatrix<Dyad,const typename M::layer_type*,bool>
edge_existence_property_matrix_undirected(
    const M* mnet
)
{
    long n = mnet->actors()->size();
    core::PropertyMatrix<Dyad,const typename M::layer_type*,bool> P(n*(n-1)/2,mnet->layers()->size(),false);

    for (auto layer: *mnet->layers())
    {
        for (auto e: *layer->edges())
        {
            Dyad d(e->v1,e->v2);
            P.set(d,layer,true);
        }
    }

    return P;
}

// @todo document: does not consider self edges
template <typename M>
core::PropertyMatrix<std::pair<const typename M::vertex_type*,const typename M::vertex_type*>,const typename M::layer_type*,bool>
edge_existence_property_matrix(
    const M* mnet
)
{
    long n = mnet->actors()->size();
    core::PropertyMatrix<std::pair<const typename M::vertex_type*,const typename M::vertex_type*>,const typename M::layer_type*,bool> P(n*(n-1),mnet->layers()->size(),false);

    for (auto layer: *mnet->layers())
    {
        if (layer->is_directed())
        {
            for (auto e: *layer->edges())
            {
                auto d = std::make_pair(e->v1,e->v2);
                P.set(d,layer,true);
            }
        }

        else
        {
            for (auto e: *layer->edges())
            {
                auto d1 = std::make_pair(e->v1,e->v2);
                auto d2 = std::make_pair(e->v2,e->v1);
                P.set(d1,layer,true);
                P.set(d2,layer,true);
            }
        }
    }

    return P;
}

// only works for multiplex networks (no inter-layer edges)
template <typename M>
core::PropertyMatrix<Triad,const typename M::layer_type*,bool>
triangle_existence_property_matrix(
    const M* mnet)
{
    long n = mnet->actors()->size();
    core::PropertyMatrix<Triad,const typename M::layer_type*,bool> P(n*(n-1)*(n-2)/6,mnet->layers()->size(),false);

    for (auto layer: *mnet->layers())
    {
        std::unordered_set<const Vertex*> processed1;

        for (auto n1: *layer->vertices())
        {
            processed1.insert(n1);
            std::unordered_set<const Vertex*> processed2;

            for (auto n2: *layer->edges()->neighbors(n1,EdgeMode::INOUT))
            {
                if (processed1.count(n2)>0)
                {
                    continue;
                }

                processed2.insert(n2);

                for (auto n3: *layer->edges()->neighbors(n2,EdgeMode::INOUT))
                {
                    if (processed1.count(n3)>0)
                    {
                        continue;
                    }

                    if (processed2.count(n3)>0)
                    {
                        continue;
                    }

                    if (layer->edges()->get(n3,n1))
                    {
                        Triad t(n1,n2,n3);
                        P.set(t,layer,true);
                    }
                }
            }
        }
    }

    return P;
}

template <typename M>
core::PropertyMatrix<const Vertex*,const typename M::layer_type*,double>
actor_degree_property_matrix(
    const M* mnet,
    EdgeMode mode
)
{
    core::PropertyMatrix<const Vertex*,const typename M::layer_type*,double> P(mnet->actors()->size(),mnet->layers()->size(),0);

    for (const Vertex* actor: *mnet->actors())
    {
        for (auto layer: *mnet->layers())
        {
            if (!layer->vertices()->contains(actor))
            {
                P.set_na(actor,layer);
            }

            else
            {
                P.set(actor,layer,layer->edges()->neighbors(actor,mode)->size());
            }
        }
    }

    return P;
}


template <typename M> core::PropertyMatrix<const Vertex*,const typename M::layer_type*,double>
actor_cc_property_matrix(
    const M* mnet
)
{
    core::PropertyMatrix<const Vertex*,const typename M::layer_type*,double> P(mnet->actors()->size(),mnet->layers()->size(),0);

    for (const Vertex* actor: *mnet->actors())
    {
        for (auto layer: *mnet->layers())
        {
            if (!layer->vertices()->contains(actor))
            {
                P.set_na(actor,layer);
            }

            else
            {
                P.set(actor,layer,cc(layer,actor));
            }
        }
    }

    return P;
}


}
}


// TODO test collisions
namespace std {
template <>
struct hash<std::pair<const uu::net::Vertex*,const uu::net::Vertex*>>
{
    size_t
    operator()(const std::pair<const uu::net::Vertex*,const uu::net::Vertex*>& d) const
    {
        size_t seed = 0;

        seed ^= hash<const uu::net::Vertex*>()(d.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<const uu::net::Vertex*>()(d.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);


        return seed;
    }
};
}

#endif
