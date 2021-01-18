#ifndef L1_HPP
#define L1_HPP

#include <RcppArmadillo.h>

namespace L1 {

    double ff(const colvec& ee) {
        return sum(abs(ee));
    }

    template <typename T> colvec ffGrd(const T& A, const colvec& s) {
        colvec grd = -A.t() * s;
        return grd;
    }

    template <typename T> int updategrd(const T& A, const colvec& s, const colvec& ee, colvec& grad, umat& ord, int j) {
        size_t i;
        bool ch = false;
        colvec u = sign(ee) - s;
        for (i = 0; i < u.size(); i++) {
            if (u[i] != 0) {
                grad -= A.row(i).t() * u[i];
                ch = true;
            }
        }
        if (ch) {
//            Rprintf("MaxAbsDeltaGrad: %f\n", max(abs(grad - ffGrd(A, sign(ee)))));
            ord = stable_sort_index(abs(grad), "descend");
            j = -1;
        }
        return j;
    }

}

#endif // L1_HPP
