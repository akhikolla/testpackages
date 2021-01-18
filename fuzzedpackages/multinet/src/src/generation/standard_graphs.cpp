#include "generation/standard_graphs.hpp"

#include "core/exceptions/assert_not_null.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include "core/utils/names.hpp"
#include "operations/_impl/add_predefined_subgraphs.hpp"
#include "objects/Vertex.hpp"

namespace uu {
namespace net {


std::unique_ptr<Network>
null_graph(
    size_t n,
    EdgeDir dir,
    bool allows_loops
)
{
    std::string name = "N_" + std::to_string(n);
    auto g = std::make_unique<Network>(name, dir, allows_loops);
    add_vertices(g.get(), n);
    return g;
}

std::unique_ptr<Network>
complete_graph(
    size_t n,
    EdgeDir dir
)
{
    std::string name = "K_" + std::to_string(n);
    bool allows_loops = false;
    auto g = std::make_unique<Network>(name, dir, allows_loops);
    add_complete_subgraph(g.get(), n);
    return g;
}


std::unique_ptr<Network>
complete_bipartite_graph(
    size_t n1,
    size_t n2,
    EdgeDir dir
)
{
    std::string name = "K_" + std::to_string(n1) + "_" + std::to_string(n2);
    bool allows_loops = false;
    auto g = std::make_unique<Network>(name, dir, allows_loops);
    add_complete_bipartite_subgraph(g.get(), n1, n2);
    return g;
}


std::unique_ptr<Network>
path_graph(
    size_t n,
    EdgeDir dir
)
{
    std::string name = "P_" + std::to_string(n);
    bool allows_loops = false;
    auto g = std::make_unique<Network>(name, dir, allows_loops);
    add_path(g.get(), n);
    return g;
}

std::unique_ptr<Network>
cycle_graph(
    size_t n,
    EdgeDir dir
)
{
    std::string name = "C_" + std::to_string(n);
    bool allows_loops = false;
    auto g = std::make_unique<Network>(name, dir, allows_loops);
    add_cycle(g.get(), n);
    return g;
}


std::unique_ptr<Network>
wheel_graph(
    size_t n
)
{
    std::string name = "W_" + std::to_string(n);
    EdgeDir dir = EdgeDir::UNDIRECTED;
    bool allows_loops = false;
    auto g = std::make_unique<Network>(name, dir, allows_loops);
    add_wheel(g.get(), n);
    return g;
}


}
}

