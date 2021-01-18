#ifndef L2_HPP
#define L2_HPP

#include <RcppArmadillo.h>

namespace L2 {

    double ff(const colvec& ee) {
         return sum(ee % ee);
     }

    template <typename T> colvec ffGrd(const T& A, const colvec& ee) {
        colvec grd = -2.0 * A.t() * ee;
        return grd;
    }

}

#endif // L2_HPP
