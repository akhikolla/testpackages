/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#include "io/_impl/read_common.hpp"
#include "core/utils/string.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include "core/exceptions/WrongFormatException.hpp"


namespace uu {
namespace net {

void
read_graph_type(
    const std::string& graph_type_spec,
    GraphType& graph_type,
    size_t line_number
)
{
    std::string feat = graph_type_spec;
    core::to_upper_case(feat);

    if (feat=="MULTI")
    {
        graph_type.allows_multi_edges=true;
    }

    else if (feat=="SIMPLE")
    {
    }

    else if (feat=="DIRECTED")
    {
        graph_type.is_directed=true;
    }

    else if (feat=="UNDIRECTED")
    {
    }

    else if (feat=="WEIGHTED")
    {
        graph_type.is_weighted=true;
    }

    else if (feat=="PROBABILISTIC")
    {
        graph_type.is_probabilistic=true;
    }


    else if (feat=="UNWEIGHTED")
    {
    }

    else if (feat=="TEMPORAL")
    {
        graph_type.is_temporal=true;
    }

    else if (feat=="STATIC")
    {
    }

    else if (feat=="LOOPS")
    {
        graph_type.allows_loops=true;
    }

    else if (feat=="NO LOOPS")
    {
    }

    else
    {
        throw core::WrongParameterException("Line " +
                                            std::to_string(line_number) +
                                            ": " + graph_type_spec);
    }
}

GraphMetadata
read_metadata(
    const std::string& infile,
    char separator
)
{

    GraphMetadata meta;

    // Set up CSV Reader

    core::CSVReader csv;
    csv.trim_fields(true);
    csv.set_field_separator(separator);
    csv.set_comment("--");

    std::string version;

    GraphIOFileSection section = GraphIOFileSection::DEFAULT;

    csv.open(infile);

    while (csv.has_next())
    {
        std::vector<std::string> fields = csv.get_next();
        std::string line = csv.get_current_raw_line();
        // remove trailing spaces
        line.erase(line.find_last_not_of(" \t")+1);
        line.erase(0,line.find_first_not_of(" \t"));

        if (line.size()==0)
        {
            continue;
        }


        if (new_section_start(line))
        {
            section = get_section(line);
            fields = csv.get_next();
            line = csv.get_current_raw_line();
            // remove trailing spaces
            line.erase(line.find_last_not_of(" \t")+1);
            line.erase(0,line.find_first_not_of(" \t"));
        }


        switch (section)
        {
        case GraphIOFileSection::VERSION:
        {
            version = read_version(line, csv.row_num());
            break;
        }

        case GraphIOFileSection::TYPE:
        {
            for (std::string graph_characteristic: fields)
            {
                read_graph_type(graph_characteristic, meta.features, csv.row_num());
            }

            break;
        }

        case GraphIOFileSection::VERTEX_ATTRIBUTES:
        {
            size_t from_idx = 0;
            core::Attribute vertex_att = read_attr_def(fields, from_idx, csv.row_num());
            meta.vertex_attributes.push_back(vertex_att);
            break;
        }

        case GraphIOFileSection::EDGE_ATTRIBUTES:
        {
            size_t from_idx = 0;
            core::Attribute edge_att = read_attr_def(fields, from_idx, csv.row_num());
            meta.edge_attributes.push_back(edge_att);
            break;
        }

        default:
            break;
        }
    }

    csv.close();
    return meta;
}




std::string
read_version(
    const std::string& line,
    size_t line_number
)
{

    std::string version = line;

    if (version != "1.0" && version != "2.0" && version != "2" && version != "3" && version != "3.0")
        throw core::WrongFormatException("Line " + std::to_string(line_number) +
                                         ": invalid version number "
                                         + line);

    // standardization
    if (version=="3")
    {
        version = "3.0";
    }

    return version;
}




core::Attribute
read_attr_def(
    const std::vector<std::string>& line,
    size_t from_idx,
    size_t line_number
)
{
    if (line.size()!=2+from_idx)
        throw core::WrongFormatException("Line " + std::to_string(line_number) +
                                         ": attribute name and type expected");

    std::string attr_name = line[from_idx+0];
    std::string attr_type_name = line[from_idx+1];
    core::to_upper_case(attr_type_name);
    core::AttributeType attr_type;

    if (attr_type_name=="NUMERIC")
    {
        attr_type = core::AttributeType::DOUBLE;    // back compatibility
    }

    else if (attr_type_name=="DOUBLE")
    {
        attr_type = core::AttributeType::DOUBLE;
    }

    else if (attr_type_name=="INT")
    {
        attr_type = core::AttributeType::INTEGER;
    }

    else if (attr_type_name=="STRING")
    {
        attr_type = core::AttributeType::STRING;
    }

    else if (attr_type_name == "TIME")
    {
        attr_type = core::AttributeType::TIME;
    }

    else throw core::WrongFormatException("Line " + std::to_string(line_number) +
                                              ": unsupported attribute type " + line[from_idx+1]);

    return core::Attribute(attr_name,attr_type);
}


}
}

