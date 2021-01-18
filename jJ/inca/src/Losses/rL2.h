#ifndef RL2_HPP
#define RL2_HPP

#include <RcppArmadillo.h>

namespace rL2 {

    double ff(const colvec& ee, const colvec& scale) {
        colvec r = ee % scale;
        return sum(r % r);
    }

    template <typename T> colvec ffGrd(const T& A, const colvec& ee, const colvec& scale) {
        colvec r = ee % scale;
        colvec grd = -2.0 * A.t() * r;
        return grd;
    }

}

#endif // RL2_HPP

