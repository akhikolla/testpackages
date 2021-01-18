#ifndef UU_MEASURES_RELEVANCE_H_
#define UU_MEASURES_RELEVANCE_H_

#include <vector>
#include "core/exceptions/assert_not_null.hpp"
#include "objects/Vertex.hpp"
#include "objects/EdgeMode.hpp"
#include "networks/_impl/containers/GenericObjectList.hpp"
#include "measures/neighborhood.hpp"

namespace uu {
namespace net {


template <typename M, typename LayerIterator>
double
relevance(
    const M* mnet,
    LayerIterator first,
    LayerIterator last,
    const Vertex* actor,
    EdgeMode mode
);

template <typename M, typename LayerIterator>
double
xrelevance(
    const M* mnet,
    LayerIterator first,
    LayerIterator last,
    const Vertex* actor,
    EdgeMode mode
);


template <typename M, typename LayerIterator>
double
relevance(
    const M* mnet,
    LayerIterator first,
    LayerIterator last,
    const Vertex* actor,
    EdgeMode mode
)
{
    double num_selected_neighbors = neighbors(first, last, actor, mode).size();
    double num_all_neighbors = neighbors(mnet->layers()->begin(), mnet->layers()->end(), actor, mode).size();

    if (num_all_neighbors==0)
    {
        return 0;    // by definition
    }

    else
    {
        return num_selected_neighbors/num_all_neighbors;
    }
}


template <typename M, typename LayerIterator>
double
xrelevance(
    const M* mnet,
    LayerIterator first,
    LayerIterator last,
    const Vertex* actor,
    EdgeMode mode
)
{
    double num_selected_neighbors = xneighbors(mnet, first, last, actor, mode).size();
    double num_all_neighbors = neighbors(mnet->layers()->begin(), mnet->layers()->end(), actor, mode).size();

    if (num_all_neighbors==0)
    {
        return 0;    // by definition
    }

    else
    {
        return num_selected_neighbors/num_all_neighbors;
    }
}

}
}

#endif
