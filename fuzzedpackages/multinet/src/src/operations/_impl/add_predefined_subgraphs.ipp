/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#include "core/exceptions/assert_not_null.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include "core/exceptions/DuplicateElementException.hpp"
#include "core/utils/names.hpp"
#include "objects/Vertex.hpp"

namespace uu {
namespace net {


template<typename G>
std::vector<const Vertex*>
add_vertices(
    G* g,
    size_t n,
    const std::string& base_vertex_name
)
{
    core::assert_not_null(g, "add_vertices", "g");

    // result vector, where to save the new vertices
    std::vector<const Vertex*> new_vertices;
    new_vertices.reserve(n);

    // adding the new vertices
    core::NameIterator names(base_vertex_name, n);

    for (auto vertex_name: names)
    {
        auto v = g->vertices()->add(vertex_name);

        if (!v)
        {
            std::string err = "Vertex " + vertex_name;
            throw core::DuplicateElementException(err);
        }

        new_vertices.push_back(v);
    }

    return new_vertices;
}


template<typename G, typename ForwardIterator>
void
add_vertices(
    G* g,
    ForwardIterator first,
    ForwardIterator last
)
{
    core::assert_not_null(g, "add_vertices", "g");

    auto current = first;

    while (current != last)
    {
        auto vertex = (*current);
        auto v = g->vertices()->add(vertex);

        if (!v)
        {
            std::string err = "Vertex " + vertex->name;
            throw core::DuplicateElementException(err);
        }

        ++current;
    }
}


template<typename G>
void
add_complete_subgraph(
    G* g,
    size_t n,
    const std::string& base_vertex_name
)
{
    auto vertices = add_vertices(g, n, base_vertex_name);

    bool dir = g->is_directed();

    for (auto v1: vertices)
    {
        for (auto v2: vertices)
        {
            if (v1 == v2)
            {
                continue;
            }

            if (!dir && v1 > v2)
            {
                continue;
            }

            g->edges()->add(v1, v2);
        }
    }
}


template<typename G>
void
add_complete_bipartite_subgraph(
    G* g,
    size_t n1,
    size_t n2,
    const std::string& base_vertex_name1,
    const std::string& base_vertex_name2
)
{
    auto vs1 = add_vertices(g, n1, base_vertex_name1);
    auto vs2 = add_vertices(g, n2, base_vertex_name2);

    bool dir = g->is_directed();

    for (auto v1: vs1)
    {
        for (auto v2: vs2)
        {
            g->edges()->add(v1, v2);

            if (dir)
            {
                g->edges()->add(v2, v1);
            }
        }
    }
}



template<typename G>
void
add_path(
    G* g,
    size_t n,
    const std::string& base_vertex_name
)
{
    auto vertices = add_vertices(g, n, base_vertex_name);

    for (size_t i = 0; i < n-1; i++)
    {
        auto v1 = vertices.at(i);
        auto v2 = vertices.at(i+1);
        g->edges()->add(v1, v2);
    }
}


template<typename G>
void
add_cycle(
    G* g,
    size_t n,
    const std::string& base_vertex_name
)
{
    auto vertices = add_vertices(g, n, base_vertex_name);

    for (size_t i = 0; i < n-1; i++)
    {
        auto v1 = vertices.at(i);
        auto v2 = vertices.at(i+1);
        g->edges()->add(v1, v2);
    }

    auto v1 = vertices.at(n-1);
    auto v2 = vertices.at(0);
    g->edges()->add(v1, v2);
}


template<typename G>
void
add_wheel(
    G* g,
    size_t n,
    const std::string& base_vertex_name
)
{
    auto vertices = add_vertices(g, n, base_vertex_name);

    auto vc = vertices.at(0);

    for (size_t i = 1; i < n-1; i++)
    {
        auto v1 = vertices.at(i);
        auto v2 = vertices.at(i+1);
        g->edges()->add(v1, v2);
        g->edges()->add(vc, v1);
    }

    auto v1 = vertices.at(n-1);
    auto v2 = vertices.at(1);
    g->edges()->add(v1, v2);
    g->edges()->add(vc, v1);
}


}
}

