#include "core/exceptions/assert_not_null.hpp"
#include "generation/empty_copy.hpp"

namespace uu {
namespace net {

template<typename G>
std::unique_ptr<G>
graph_union(
    const G* g1,
    const G* g2,
    const std::string& name
)
{
    core::assert_not_null(g1, "graph_union", "g1");
    core::assert_not_null(g2, "graph_union", "g2");

    if (g1->is_directed() != g2->is_directed())
    {
        std::string err = "union between directed and undirected graphs not supported";
        throw core::OperationNotSupportedException(err);
    }

    std::unique_ptr<G> res = empty_copy(g1, name);

    // computing the union of the vertices and edges

    std::set<const Vertex*> v_union;

    for (auto vertex: *g1->vertices())
    {
        v_union.insert(vertex);
    }

    for (auto vertex: *g2->vertices())
    {
        v_union.insert(vertex);
    }

    std::set<const Edge*> e_union;

    for (auto edge: *g1->edges())
    {
        e_union.insert(edge);
    }

    for (auto edge: *g2->edges())
    {
        e_union.insert(edge);
    }

    // adding the unions of vertices and edges to the new graph

    for (auto vertex: v_union)
    {
        auto success = res->vertices()->add(vertex);

        if (!success)
        {
            std::string err = "the two networks contain different vertices with the same name";
            throw core::OperationNotSupportedException(err);
        }
    }

    for (auto edge: e_union)
    {
        auto success = res->edges()->add(edge);

        if (!success)
        {
            std::string err = "the two networks contain different edges with the same end-vertices";
            throw core::OperationNotSupportedException(err);
        }
    }

    return res;
}

template<typename G>
void
graph_add(
    const G* g,
    G* target)
{

    core::assert_not_null(g, "graph_add", "g");
    core::assert_not_null(target, "graph_add", "target");


    for (auto vertex: *g->vertices())
    {
        target->vertices()->add(vertex);
    }


    for (auto edge: *g->edges())
    {
        target->edges()->add(edge->v1,edge->v2);
    }

    if (!g->is_directed() && target->is_directed())
    {
        for (auto edge: *g->edges())
        {
            target->edges()->add(edge->v2,edge->v1);
        }
    }
}

template<typename G, typename W>
void
weigthed_graph_add(
    const G* g,
    W* target)
{

    core::assert_not_null(g, "weigthed_graph_add", "g");
    core::assert_not_null(target, "weigthed_graph_add", "target");


    for (auto vertex: *g->vertices())
    {
        target->vertices()->add(vertex);
    }


    for (auto edge: *g->edges())
    {
        auto new_edge = target->edges()->get(edge->v1,edge->v2);

        if (!new_edge)
        {
            new_edge = target->edges()->add(edge->v1,edge->v2);
            target->edges()->attr()->set_weight(new_edge,1.0);
        }

        else
        {
            double weight = target->edges()->attr()->get_weight(new_edge).value;
            target->edges()->attr()->set_weight(new_edge, weight+1.0);
        }
    }

    if (!g->is_directed() && target->is_directed())
    {
        for (auto edge: *g->edges())
        {
            auto new_edge = target->edges()->get(edge->v2,edge->v1);

            if (!new_edge)
            {
                new_edge = target->edges()->add(edge->v2,edge->v1);
                target->edges()->attr()->set_weight(new_edge,1.0);
            }

            else
            {
                double weight = target->edges()->attr()->get_weight(new_edge).value;
                target->edges()->attr()->set_weight(new_edge, weight+1.0);
            }
        }
    }

}



template<typename GraphIterator, typename G>
void
weighted_graph_union(
    GraphIterator begin,
    GraphIterator end,
    G* target,
    const std::string& weight_attribute_name
)
{
    core::assert_not_null(target, "graph_union", "target");

    for (auto g=begin; g!=end; ++g)
    {
        weigthed_graph_add(*g, target, weight_attribute_name);
    }
}

template<typename G, typename W>
void
weigthed_graph_add(
    const G* g,
    W* target,
    const std::string& weight_attribute_name
)
{

    core::assert_not_null(g, "weigthed_graph_add", "g");
    core::assert_not_null(target, "weigthed_graph_add", "target");


    for (auto vertex: *g->vertices())
    {
        target->vertices()->add(vertex);
    }


    for (auto edge: *g->edges())
    {
        auto new_edge = target->edges()->get(edge->v1,edge->v2);

        if (!new_edge)
        {
            new_edge = target->edges()->add(edge->v1,edge->v2);
            target->edges()->attr()->set_double(new_edge, weight_attribute_name,1.0);
        }

        else
        {
            double weight = target->edges()->attr()->get_double(new_edge, weight_attribute_name).value;
            target->edges()->attr()->set_double(new_edge, weight_attribute_name, weight+1.0);
        }
    }

    if (!g->is_directed() && target->is_directed())
    {
        for (auto edge: *g->edges())
        {
            auto new_edge = target->edges()->get(edge->v2,edge->v1);

            if (!new_edge)
            {
                new_edge = target->edges()->add(edge->v2,edge->v1);
                target->edges()->attr()->set_double(new_edge, weight_attribute_name, 1.0);
            }

            else
            {
                double weight = target->edges()->attr()->get_double(new_edge, weight_attribute_name).value;
                target->edges()->attr()->set_double(new_edge, weight_attribute_name, weight+1.0);
            }
        }
    }

}



}
}

