#ifndef AL1_HPP
#define AL1_HPP

#include <RcppArmadillo.h>

namespace aL1 {

    double ff(const colvec& ee, const colvec& lambda) {
        return sum(ee % lambda); // here "colvec % colvec" is a dot product (while "int % int" is the <mod> operator)
    }

    template <typename T> colvec ffGrd(const T& A, const colvec& lambda) {
        colvec grd = -A.t() * lambda;
        return grd;
    }

    template <typename T> int updategrd(const T& A, const colvec& s, const colvec& ee, colvec& grad, umat& ord, int j) {
        size_t i;
        int u;
        if (any((s < 0) != (ee < 0))) {
            for (i = 0; i < s.size(); i++) {
                u = (ee[i] < 0) - (s[i] < 0);
                if (u != 0) {
                    grad -= A.row(i).t() * u;
                }
            }
            ord = stable_sort_index(abs(grad), "descend");
            j = -1;
        }
        return j;
    }

}

#endif // AL1_HPP
