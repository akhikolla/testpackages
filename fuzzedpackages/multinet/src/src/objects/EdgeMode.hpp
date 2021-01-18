#ifndef UU_OBJECTS_EDGEMODE_H_
#define UU_OBJECTS_EDGEMODE_H_

namespace uu {
namespace net {

/**
 * Selection mode, for directed edges.
 * E.g., to compute the IN-degree or OUT-degree of a node.
 */
enum class EdgeMode
{
    INOUT,
    IN,
    OUT
};

}
}

#endif
