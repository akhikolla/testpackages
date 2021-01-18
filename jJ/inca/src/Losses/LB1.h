#ifndef LB1_HPP
#define LB1_HPP

#include <RcppArmadillo.h>

namespace LB1 {

    double ff(const colvec& L, const colvec& U, const colvec& e) {
        size_t i;
        double loss = 0.0;
        for (i = 0; i < e.size(); i++) {
            if (L[i] > e[i]) {
                loss += fabs(e[i] - L[i]);
            }
            if (U[i] < e[i]) {
                loss += fabs(U[i] - e[i]);
            }
        }
        return loss;
    }

    template <typename T> colvec ffGrd(const T& A, const colvec& e, const colvec& L, const colvec& U) {
        colvec Df = conv_to<colvec>::from(e > U) - conv_to<colvec>::from(e < L);
        colvec grd = - A.t() * Df;
        return grd;
    }

    template <typename T> int updategrd(const T& A, const mat& B, const colvec& s, const colvec& ee, colvec& grad, umat& ord, int j) {
        size_t i;
        bool ch = false;
        colvec u = conv_to<colvec>::from(ee > B.col(1)) - conv_to<colvec>::from(ee < B.col(0));
        u -= s;
        for (i = 0; i < u.size(); i++) {
            if (u[i] != 0) {
                grad -= A.row(i).t() * u[i];
                ch = true;
            }
        }
        if (ch) {
            ord = stable_sort_index(abs(grad), "descend");
            j = -1;
        }
        return j;

    }

}

#endif // LB1_HPP

