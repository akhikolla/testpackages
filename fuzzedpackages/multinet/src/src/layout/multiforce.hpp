#ifndef UU_LAYOUT_MULTIFORCE_H_
#define UU_LAYOUT_MULTIFORCE_H_

#include <map>
#include <unordered_map>
#include "core/utils/random.hpp"
#include "core/exceptions/assert_not_null.hpp"
#include "layout/XYZCoordinates.hpp"
#include "objects/Vertex.hpp"

namespace uu {
namespace net {

double
fr(double p, double k);

double
fain(double p, double k);

double
fainter(double p, double k);

template <typename M>
std::map<std::pair<const Vertex*, const typename M::layer_type*>,XYZCoordinates>
multiforce(
    const M* mnet,
    double width,
    double length,
    const std::unordered_map<const typename M::layer_type*,double>& w_intra,
    const std::unordered_map<const typename M::layer_type*,double>& w_inter,
    const std::unordered_map<const typename M::layer_type*,double>& gravity,
    int iterations
);



}
}

#include "multiforce.ipp"

#endif
