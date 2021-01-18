#ifndef UU_OBJECTS_COMPARISONRESULT_H_
#define UU_OBJECTS_COMPARISONRESULT_H_

namespace uu {
namespace net {

/**
 * Outcome of the comparison between two numbers, or two vectors, or two matrices.
 * (the concept of Pareto dominance is used).
 */
enum class ComparisonResult
{
    LESS_THAN,
    EQUAL,
    INCOMPARABLE,
    GREATER_THAN
};


}
}

#endif
