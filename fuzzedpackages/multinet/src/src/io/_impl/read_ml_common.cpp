#include "io/_impl/read_ml_common.hpp"
#include "io/_impl/read_common.hpp"
#include "core/utils/string.hpp"
#include <iostream>
#include "core/exceptions/WrongFormatException.hpp"

namespace uu {
namespace net {

bool
new_multilayer_section_start(const std::string& line)
{
    if (!(line.find("#")==0))
    {
        return false;
    }

    std::string line_copy = line;
    core::to_upper_case(line_copy);

    if (
        line_copy.find("#VERSION") == 0 || // for backward compatibility, if the version is on the same line
        line_copy.find("#TYPE") == 0 || // for backward compatibility, if the type is on the same line
        line_copy=="#LAYERS" ||
        line_copy=="#ACTORS" ||
        line_copy=="#ACTOR ATTRIBUTES" ||
        line_copy=="#VERTICES" ||
        line_copy=="#VERTEX ATTRIBUTES" ||
        // OLD line_copy=="#INTRALAYER EDGES" ||
        // OLD line_copy=="#INTERLAYER EDGES" ||
        line_copy=="#EDGES" ||
        line_copy=="#EDGE ATTRIBUTES" ||
        // deprecated
        line_copy=="#VERTEXES")
    {
        return true;
    }

    return false;
}

MultilayerIOFileSection
get_multilayer_section(
    const std::string& line
)
{
    std::string line_copy = line;
    core::to_upper_case(line_copy);

    if (line_copy.find("#VERSION") == 0) // for backward compatibility, if the type is on the same line
    {
        return MultilayerIOFileSection::VERSION;
    }

    if (line_copy.find("#TYPE") == 0) // for backward compatibility, if the type is on the same line
    {
        return MultilayerIOFileSection::TYPE;
    }

    if (line_copy=="#LAYERS")
    {
        return MultilayerIOFileSection::LAYERS;
    }

    if (line_copy=="#ACTORS")
    {
        return MultilayerIOFileSection::ACTORS;
    }

    if (line_copy=="#ACTOR ATTRIBUTES")
    {
        return MultilayerIOFileSection::ACTOR_ATTRIBUTES;
    }

    if (line_copy=="#VERTICES" || line_copy=="#NODES")
    {
        return MultilayerIOFileSection::INTRALAYER_VERTICES;
    }

    if (line_copy=="#VERTEX ATTRIBUTES" || line_copy=="#NODE ATTRIBUTES")
    {
        return MultilayerIOFileSection::VERTEX_ATTRIBUTES;
    }

    if (line_copy=="#EDGES")
    {
        return MultilayerIOFileSection::EDGES;
    }

    if (line_copy=="#EDGE ATTRIBUTES")
    {
        return MultilayerIOFileSection::EDGE_ATTRIBUTES;
    }

    // DEPRECATED
    if (line_copy=="#VERTEXES")
    {
        // std::cerr << "[WARNING] usage of #VERTEXES deprecated. Use #VERTICES instead." << std::endl;
        return MultilayerIOFileSection::INTRALAYER_VERTICES;
    }

    return MultilayerIOFileSection::DEFAULT; // cannot get here
}


MultilayerMetadata
read_multilayer_metadata(
    const std::string& infile,
    char separator
)
{

    MultilayerMetadata meta;

    // Set up CSV Reader

    core::CSVReader csv;
    csv.trim_fields(true);
    csv.set_field_separator(separator);
    csv.set_comment("--");

    std::string version;

    MultilayerIOFileSection section = MultilayerIOFileSection::DEFAULT;

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


        if (new_multilayer_section_start(line))
        {
            section = get_multilayer_section(line);

            if (section == MultilayerIOFileSection::TYPE) // for backward compatibility
            {
                std::string line_copy = line;
                core::to_upper_case(line_copy);

                if (line_copy.find("MULTIPLEX") != std::string::npos)
                {
                    meta.is_multiplex = true;
                }

                else if (line_copy.find("MULTILAYER") != std::string::npos)
                {
                    meta.is_multiplex = false;
                }
            }

            //fields = csv.get_next();
            //line = csv.get_current_raw_line();
            // remove trailing spaces
            //line.erase(line.find_last_not_of(" \t")+1);
            //line.erase(0,line.find_first_not_of(" \t"));
            continue;
        }


        switch (section)
        {

        case MultilayerIOFileSection::TYPE:
        {
            std::string line_copy = line;
            core::to_upper_case(line_copy);

            if (line_copy == "MULTIPLEX")
            {
                meta.is_multiplex = true;
            }

            else if (line_copy == "MULTILAYER")
            {
                meta.is_multiplex = false;
            }

            else
            {
                throw core::WrongFormatException("Line " + std::to_string(csv.row_num()) +
                                                 ": wrong network type (" + line + ")");
            }

            break;
        }

        case MultilayerIOFileSection::VERSION:
        {
            version = read_version(line, csv.row_num());

            if (version != "3.0")
                throw core::WrongFormatException("Line " + std::to_string(csv.row_num()) +
                                                 ": version 3.0 required");

            break;
        }


        case MultilayerIOFileSection::LAYERS:
        {
            if (meta.is_multiplex)
            {
                if (fields.size() < 2)
                {
                    throw core::WrongFormatException("Line " + std::to_string(csv.row_num()) +
                                                     ": layer name and directionality required");
                }

                std::string layer_name = fields.at(0);
                meta.layers[layer_name];

                for (size_t idx = 1; idx<fields.size(); idx++)
                {
                    read_graph_type(fields.at(idx), meta.layers[layer_name], csv.row_num());
                }
            }

            else
            {
                if (fields.size() < 3)
                {
                    throw core::WrongFormatException("Line " + std::to_string(csv.row_num()) +
                                                     ": pair of layer names and directionality required");
                }

                std::string layer_name1 = fields.at(0);
                meta.layers[layer_name1];
                std::string layer_name2 = fields.at(1);
                meta.layers[layer_name2];

                if (layer_name1 == layer_name2)
                {
                    meta.layers[layer_name1];

                    for (size_t idx = 2; idx<fields.size(); idx++)
                    {
                        read_graph_type(fields.at(idx), meta.layers[layer_name1], csv.row_num());
                    }
                }

                else
                {
                    std::string dir = fields.at(2);
                    core::to_upper_case(dir);

                    if (dir=="DIRECTED")
                    {
                        meta.interlayer_dir[std::pair<std::string,std::string>(layer_name1,layer_name2)] = true;
                    }
                }

            }

            break;
        }


        case MultilayerIOFileSection::ACTOR_ATTRIBUTES:
        {
            if (fields.size()==2)
            {
                size_t from_idx = 0;
                core::Attribute vertex_att = read_attr_def(fields, from_idx, csv.row_num());
                meta.vertex_attributes.push_back(vertex_att);
            }

            else
            {
                throw core::WrongFormatException("Line " + std::to_string(csv.row_num()) +
                                                 ": attribute name and attribute type expected");
            }

            break;
        }


        case MultilayerIOFileSection::VERTEX_ATTRIBUTES:
        {
            // intralayer vertex attributes
            if (fields.size()==3)
            {
                // add layer if not previously defined
                std::string layer_name = fields[0];

                if (meta.layers.find(layer_name) == meta.layers.end())
                {
                    GraphType gt;
                    meta.layers[layer_name] = gt;
                }

                // read attribute
                size_t from_idx = 1;
                core::Attribute vertex_att = read_attr_def(fields, from_idx, csv.row_num());
                meta.intralayer_vertex_attributes[layer_name].push_back(vertex_att);
            }

            else
            {
                throw core::WrongFormatException("Line " + std::to_string(csv.row_num()) +
                                                 ": layer, attribute name and attribute type expected");
            }

            break;
        }

        case MultilayerIOFileSection::EDGE_ATTRIBUTES:
        {
            if (fields.size()==2)
            {
                int from_idx = 0;
                core::Attribute edge_att = read_attr_def(fields, from_idx, csv.row_num());
                meta.interlayer_edge_attributes.push_back(edge_att);
            }

            else if (fields.size()==3)
            {
                // add layer if not previously defined
                std::string layer_name = fields[0];

                if (meta.layers.find(layer_name) == meta.layers.end())
                {
                    GraphType gt;
                    meta.layers[layer_name] = gt;
                }

                // read attribute
                int from_idx = 1;
                core::Attribute edge_att = read_attr_def(fields, from_idx, csv.row_num());
                meta.intralayer_edge_attributes[layer_name].push_back(edge_att);
            }

            else
            {
                throw core::WrongFormatException("Line " + std::to_string(csv.row_num()) +
                                                 ": layer name (optional), attribute name and attribute type expected");
            }

            break;
        }

        default:
            break;
        }
    }

    // @todo global edge attributes added to all layers at the end. Maybe add more flexibility
    for (auto edge_att: meta.interlayer_edge_attributes)
    {
        for (auto layer: meta.layers)
        {
            std::string layer_name = layer.first;
            meta.intralayer_edge_attributes[layer_name].push_back(edge_att);
        }
    }

    csv.close();
    return meta;
}

}
}
