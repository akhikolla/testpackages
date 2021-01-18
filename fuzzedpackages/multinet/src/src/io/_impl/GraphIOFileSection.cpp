/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#include <iostream>
#include "io/_impl/GraphIOFileSection.hpp"
#include "core/utils/string.hpp"
#include "core/exceptions/WrongParameterException.hpp"


namespace uu {
namespace net {

bool
new_section_start(const std::string& line)
{
    if (!(line.find("#")==0))
    {
        return false;
    }

    std::string line_copy = line;
    core::to_upper_case(line_copy);

    if (
        line_copy=="#VERSION" ||
        line_copy=="#TYPE" ||
        line_copy=="#VERTEX ATTRIBUTES" ||
        line_copy=="#EDGE ATTRIBUTES" ||
        line_copy=="#VERTICES" ||
        line_copy=="#EDGES" ||
        // deprecated
        line_copy=="#VERTEXES" ||
        line_copy=="#ACTORS" ||
        line_copy=="#ACTOR ATTRIBUTES")

    {
        return true;
    }

    return false;
}



GraphIOFileSection
get_section(const std::string& line)
{
    std::string line_copy = line;
    core::to_upper_case(line_copy);

    if (line_copy=="#VERSION")
    {
        return GraphIOFileSection::VERSION;
    }

    if (line_copy=="#TYPE")
    {
        return GraphIOFileSection::TYPE;
    }

    if (line_copy=="#VERTEX ATTRIBUTES")
    {
        return GraphIOFileSection::VERTEX_ATTRIBUTES;
    }

    if (line_copy=="#EDGE ATTRIBUTES")
    {
        return GraphIOFileSection::EDGE_ATTRIBUTES;
    }

    if (line_copy=="#VERTICES")
    {
        return GraphIOFileSection::VERTICES;
    }

    if (line_copy=="#EDGES")
    {
        return GraphIOFileSection::EDGES;
    }

    // DEPRECATED
    if (line_copy=="#VERTEXES")
    {
        //std::cerr << "[WARNING] usage of #VERTEXES deprecated. Use #VERTICES instead." << std::endl;
        return GraphIOFileSection::VERTICES;
    }

    if (line_copy=="#ACTORS")
    {
        //std::cerr << "[WARNING] usage of #ACTORS deprecated. Use #VERTICES instead." << std::endl;
        return GraphIOFileSection::VERTICES;
    }

    if (line_copy=="#ACTOR ATTRIBUTES")
    {
        //std::cerr << "[WARNING] usage of #ACTOR deprecated. Use #VERTEX instead." << std::endl;
        return GraphIOFileSection::VERTEX_ATTRIBUTES;
    }

    return GraphIOFileSection::DEFAULT; // cannot get here
}


}
}

