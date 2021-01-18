#include "io/read_network.hpp"

namespace uu {
namespace net {

std::unique_ptr<Network>
read_network(
    const std::string& infile,
    const std::string& name,
    char separator
)
{
    // Read metadata
    GraphMetadata meta = read_metadata(infile, ',');
    EdgeDir dir = meta.features.is_directed?EdgeDir::DIRECTED:EdgeDir::UNDIRECTED;

    // Check metadata consistency (@todo)
    // create network
    // and add attributes
    auto g = std::make_unique<Network>(name, dir, meta.features.allows_loops);

    for (auto attr: meta.vertex_attributes)
    {
        g->vertices()->attr()->add(attr.name, attr.type);
    }

    for (auto attr: meta.edge_attributes)
    {
        g->edges()->attr()->add(attr.name, attr.type);
    }

    // Read data (vertices, edges, attribute values)
    read_data(g.get(),  meta, infile, separator);

    return g;

}


template <>
void
read_vertex(
    Network* g,
    const std::vector<std::string>& fields,
    const std::vector<core::Attribute>& vertex_attributes,
    size_t line_number
)
{


    if (fields.size()>1+vertex_attributes.size())
    {
        throw core::WrongFormatException("Line " +
                                         std::to_string(line_number) +
                                         ": wrong number of attributes");
    }

    auto v = read_vertex(g, fields, 0, line_number);


    read_attr_values(g->vertices()->attr(), v, vertex_attributes, fields, 1, line_number);


}


template <>
void
read_edge(
    Network* g,
    const std::vector<std::string>& fields,
    const std::vector<core::Attribute>& edge_attributes,
    size_t line_number
)
{

    if (fields.size()!=2+edge_attributes.size())
    {
        throw core::WrongFormatException("Line " +
                                         std::to_string(line_number) +
                                         ": From and To actor names and weight must " +
                                         "be specified for each edge");
    }

    auto edge = read_edge(g, fields, 0, line_number);

    read_attr_values(g->edges()->attr(), edge, edge_attributes, fields, 2, line_number);

}

}
}
