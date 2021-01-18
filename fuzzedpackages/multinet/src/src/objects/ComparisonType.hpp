#ifndef UU_OBJECTS_COMPARISONTYPE_H_
#define UU_OBJECTS_COMPARISONTYPE_H_

namespace uu {
namespace net {

/**
 * The length of a path traversing multiple layers can be represented
 * using a matrix m, where cell m_{i,j} indicates the number
 * of edges traversed from layer i to layer j.
 * These constants specify how much of this information should be used
 * while comparing two path lengths:
 * - FULL: all the elements in the matrix are considered.
 * - SWITCH_COSTS: the diagonal elements are considered (indicating inter-layer steps),
 *   plus an additional element summing the values on all non-diagonal cells (indicating a layer crossing).
 * - MULTIPLEX: only the diagonal is considered (that is, only inter-layer steps).
 * - SIMPLE: the sum of all elements in the matrix (that is, the number of steps
 *   irrespective of the traversed layers) is considered.
 **/
enum class ComparisonType
{
    FULL,
    SWITCH_COSTS,
    MULTIPLEX,
    SIMPLE
};



}
}

#endif
