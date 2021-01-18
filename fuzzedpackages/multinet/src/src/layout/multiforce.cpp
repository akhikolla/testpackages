
/**
 * History:
 * - 2018.06.22 replaced rescaling of the node positions with hard frames, as in the F/R algorithm.
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */
#include "layout/multiforce.hpp"

namespace uu {
namespace net {

/** repulsive force */
double
fr(double d, double k)
{
    return k*k/d;
}

/** attractive force, inter-layer */
double
fain(double d, double k)
{
    return d*d/k;
}

/** attractive force, intra-layer */
double
fainter(double d, double k)
{
    return d*d/k;
}

}
}


