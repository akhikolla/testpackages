#ifndef RL1_HPP
#define RL1_HPP

#include <RcppArmadillo.h>

namespace rL1 {

    double ff(const colvec& ee, const colvec& scale) {
        return sum(abs(ee) % scale);
    }

    template <typename T> colvec ffGrd(const T& A, const colvec& s, const colvec& scale) {
        colvec grd = -A.t() * (s % scale);
        return grd;
    }

    template <typename T> int updategrd(const T& A, const colvec& scale, const colvec& s, const colvec& ee, colvec& grad, umat& ord, int j) {
        size_t i;
        colvec u = sign(ee) - s;
        if (any(u != 0)) {
            for (i = 0; i < u.size(); i++) {
                if (u[i] != 0) {
                    grad -= A.row(i).t() * (u[i] * scale[i]);
                }
            }
            ord = stable_sort_index(abs(grad), "descend");
            j = -1;
        }
        return j;
    }

}

#endif // RL1_HPP

