#ifndef UU_MEASURES_MLDEGREE_H_
#define UU_MEASURES_MLDEGREE_H_

#include <vector>
#include "core/exceptions/assert_not_null.hpp"
#include "core/utils/math.hpp"
#include "objects/Vertex.hpp"
#include "objects/EdgeMode.hpp"
#include "measures/degree.hpp"

namespace uu {
namespace net {

/**
 * Returns the sum of the intralayer degrees of a vertex.
 * @param first, last iterators specifying a range of layers (first included, last not included)
 * @param v input vertex
 * @param mode to select IN, OUT, or INOUT degree
 * @return the sum of the degrees of v in the input layers
 */
template <typename LayerIterator>
int
degree(
    LayerIterator first,
    LayerIterator last,
    const Vertex* v,
    EdgeMode mode
);

/**
 * Returns the average intralayer degree of a vertex.
 * @param first, last iterators specifying a range of layers (first included, last not included)
 * @param v input vertex
 * @param mode to select IN, OUT, or INOUT degree
 * @return the mean of the degrees of v in the input layers
 */
template <typename LayerIterator>
double
degree_mean(
    LayerIterator first,
    LayerIterator last,
    const Vertex* v,
    EdgeMode mode
);

/**
 * Returns the standard deviation of the intralayer degrees of a vertex.
 * @param first, last iterators specifying a range of layers (first included, last not included)
 * @param v input vertex
 * @param mode to select IN, OUT, or INOUT degree
 * @return the standard deviation of the degrees of v in the input layers
 */
template <typename LayerIterator>
double
degree_deviation(
    LayerIterator first,
    LayerIterator last,
    const Vertex* v,
    EdgeMode mode
);

// DEFINITIONS


template <typename LayerIterator>
int
degree(
    LayerIterator first,
    LayerIterator last,
    const Vertex* v,
    EdgeMode mode
)
{
    core::assert_not_null(v, "degree", "v");

    int d = 0;

    for (auto layer=first; layer!=last; ++layer)
    {
        d += degree(*layer, v, mode);
    }

    return d;
}


template <typename LayerIterator>
double
degree_mean(
    LayerIterator first,
    LayerIterator last,
    const Vertex* v,
    EdgeMode mode
)
{
    core::assert_not_null(v, "degree_mean", "v");

    std::vector<double> degrees;

    for (auto layer=first; layer!=last; ++layer)
    {
        degrees.push_back((double)degree(*layer, v, mode));
    }

    return core::mean(degrees.begin(), degrees.end());
}


template <typename LayerIterator>
double
degree_deviation(
    LayerIterator first,
    LayerIterator last,
    const Vertex* v,
    EdgeMode mode
)
{
    core::assert_not_null(v, "degree_deviation", "v");

    std::vector<double> degrees;

    for (auto layer=first; layer!=last; ++layer)
    {
        degrees.push_back((double)degree(*layer, v, mode));
    }

    return core::stdev(degrees.begin(), degrees.end());
}


}
}

#endif
