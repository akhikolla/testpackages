#ifndef AL2_HPP
#define AL2_HPP

#include <RcppArmadillo.h>

namespace aL2 {

    double ff(const colvec& ee, const colvec& lambda) {
        return sum(ee % ee % lambda); // here "colvec % colvec" is a dot product (while "int % int" is the <mod> operator)
    }

    template <typename T> colvec ffGrd(const T& A, const colvec& ee, const colvec& lambda) {
        colvec grd = -2.0 * A.t() * (ee % lambda);
        return grd;
    }

}

#endif // AL2_HPP

