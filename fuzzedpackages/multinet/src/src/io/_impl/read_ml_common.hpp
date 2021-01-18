/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_MNET_IO_READCOMMON_H_
#define UU_MNET_IO_READCOMMON_H_

#include <string>
#include <vector>
#include "core/exceptions/OperationNotSupportedException.hpp"
#include "core/attributes/Attribute.hpp"
//#include "core/attributes/AttributeValueMap.hpp"
#include "core/utils/CSVReader.hpp"
#include "core/datastructures/objects/Object.hpp"
#include "networks/_impl/Graph.hpp"
#include "objects/Vertex.hpp"
#include "objects/Edge.hpp"
#include "io/_impl/MultilayerMetadata.hpp"
#include "io/_impl/MultilayerIOFileSection.hpp"

namespace uu {
namespace net {

/**
 * Checks if the input line indicates the start of a new section.
 */
bool
new_multilayer_section_start(
    const std::string& line
);

/**
 * Returns the new section starting on this input line.
 */
MultilayerIOFileSection
get_multilayer_section(
    const std::string& line
);


MultilayerMetadata
read_multilayer_metadata(
    const std::string& infile,
    char separator
);


template <typename ML>
void
read_multilayer_data(
    ML* ml,
    const MultilayerMetadata& meta,
    const std::string& infile,
    char separator
);

/**
 * Utility function to read ...
 * @param graph_type...
 * @param line_number current line in the input file, for error management
 */
void
read_multilayer_type(
    const std::string& graph_type,
    MultilayerMetadata& meta,
    size_t line_number
);


template <typename ML>
const Vertex*
read_actor(
    ML* g,
    const std::vector<std::string>& fields,
    size_t from_idx,
    size_t line_number
)
{
    (void)line_number; // attribute not used

    core::assert_not_null(g, "read_actor", "g");

    std::string actor_name = fields.at(from_idx);

    auto actor = g->actors()->add(actor_name);

    if (!actor)
    {
        actor = g->actors()->get(actor_name);
    }

    return actor;
}

/*
template <typename G>
void
read_edge(
    G* g,
    const std::vector<std::string>& fields,
    const std::vector<core::Attribute>& edge_attributes,
    size_t line_number
);
*/

template <typename ML>
void
read_multilayer_data(
    ML* ml,
    const MultilayerMetadata& meta,
    const std::string& infile,
    char separator
)
{

    // PASS 2: read the data

    core::CSVReader csv;
    csv.trim_fields(true);
    csv.set_field_separator(separator);
    csv.set_comment("--");
    csv.open(infile);

    MultilayerIOFileSection section = MultilayerIOFileSection::EDGES;

    while (csv.has_next())
    {
        std::vector<std::string> fields = csv.get_next();
        std::string line = csv.get_current_raw_line();
        // std::cout << line << std::endl;
        // remove trailing spaces
        line.erase(line.find_last_not_of(" \t")+1);
        line.erase(0,line.find_first_not_of(" \t"));

        if (line.size()==0)
        {
            continue;
        }


        // if new section starts here, set the current section and proceed
        if (new_multilayer_section_start(line))
        {
            section = get_multilayer_section(line);
            //fields = csv.get_next();
            continue;
        }

        switch (section)
        {
        case MultilayerIOFileSection::ACTORS:
        {

            read_vertex(ml, fields, meta, csv.row_num());
            break;
        }

        case MultilayerIOFileSection::INTRALAYER_VERTICES:
        {
            read_intralayer_vertex(ml, fields, meta, csv.row_num());
            break;
        }

        case MultilayerIOFileSection::EDGES:
        {

            if (meta.is_multiplex)
            {
                read_intralayer_edge(ml, fields, meta, csv.row_num());
            }

            else
            {
                read_interlayer_edge(ml, fields, meta, csv.row_num());
            }

            break;
        }

        default:
            break;
        }
    }

}



template <typename ML, typename G>
G*
read_layer(
    ML* ml,
    const std::vector<std::string>& fields,
    size_t from_idx,
    size_t line_number
)
{
    (void)ml; // attribute not used
    (void)fields; // attribute not used
    (void)from_idx; // attribute not used
    (void)line_number; // attribute not used
    throw core::OperationNotSupportedException("Network type not supported (IO)");
}

template <typename ML, typename G>
G*
read_layer(
    ML* ml,
    const std::vector<std::string>& fields,
    const MultilayerMetadata& meta,
    size_t line_number
)
{
    (void)ml; // attribute not used
    (void)fields; // attribute not used
    (void)meta; // attribute not used
    (void)line_number; // attribute not used
    throw core::OperationNotSupportedException("Network type not supported (IO)");
}

template <typename ML>
void
read_vertex(
    ML* ml,
    const std::vector<std::string>& fields,
    const MultilayerMetadata& meta,
    size_t line_number
)
{
    (void)ml; // attribute not used
    (void)fields; // attribute not used
    (void)meta; // attribute not used
    (void)line_number; // attribute not used
    throw core::OperationNotSupportedException("Network type not supported (IO)");
}

template <typename ML>
void
read_intralayer_vertex(
    ML* ml,
    const std::vector<std::string>& fields,
    const MultilayerMetadata& meta,
    size_t line_number
)
{
    (void)ml; // attribute not used
    (void)fields; // attribute not used
    (void)meta; // attribute not used
    (void)line_number; // attribute not used
    throw core::OperationNotSupportedException("Network type not supported (IO)");
}

template <typename ML>
void
read_intralayer_edge(
    ML* ml,
    const std::vector<std::string>& fields,
    const MultilayerMetadata& meta,
    size_t line_number
)
{
    (void)ml; // attribute not used
    (void)fields; // attribute not used
    (void)meta; // attribute not used
    (void)line_number; // attribute not used
    throw core::OperationNotSupportedException("Network type not supported (IO)");
}

template <typename ML>
void
read_interlayer_edge(
    ML* ml,
    const std::vector<std::string>& fields,
    const MultilayerMetadata& meta,
    size_t line_number
)
{
    (void)ml; // attribute not used
    (void)fields; // attribute not used
    (void)meta; // attribute not used
    (void)line_number; // attribute not used
    throw core::OperationNotSupportedException("Network type not supported (IO)");
}


}
}

#endif
