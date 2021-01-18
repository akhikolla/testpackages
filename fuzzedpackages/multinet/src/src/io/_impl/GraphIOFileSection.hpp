/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_IO_GRAPHIOFILESECTION_H_
#define UU_NET_IO_GRAPHIOFILESECTION_H_

#include <string>

namespace uu {
namespace net {

/**
 * The sections that can be found in a file storing a graph.
 */
enum class GraphIOFileSection
{
    DEFAULT, // @todo is this necessary?
    VERSION,
    TYPE,
    VERTEX_ATTRIBUTES,
    EDGE_ATTRIBUTES,
    VERTICES,
    EDGES
};

/**
 * Checks if the input line indicates the start of a new section.
 */
bool
new_section_start(
    const std::string& line
);

/**
 * Returns the new section starting on this input line.
 */
GraphIOFileSection
get_section(
    const std::string& line
);

}
}

#endif
