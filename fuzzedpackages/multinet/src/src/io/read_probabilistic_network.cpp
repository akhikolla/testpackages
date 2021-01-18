#include "io/read_probabilistic_network.hpp"

namespace uu {
namespace net {

std::unique_ptr<ProbabilisticNetwork>
read_probabilistic_network(
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
    auto g = std::make_unique<ProbabilisticNetwork>(name, dir, meta.features.allows_loops);


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
    ProbabilisticNetwork* g,
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
    ProbabilisticNetwork* g,
    const std::vector<std::string>& fields,
    const std::vector<core::Attribute>& edge_attributes,
    size_t line_number
)
{

    if (fields.size()!=3+edge_attributes.size())
    {
        throw core::WrongFormatException("Line " +
                                         std::to_string(line_number) +
                                         ": From and To vertex names and probability must " +
                                         "be specified for each edge");
    }

    auto edge = read_edge(g, fields, 0, line_number);

    double p = core::to_double(fields.at(2));

    g->set_prob(edge, p);

    read_attr_values(g->edges()->attr(), edge, edge_attributes, fields, 3, line_number);

}

}
}
