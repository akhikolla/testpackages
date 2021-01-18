#ifndef RB1_HPP
#define RB1_HPP

#include <RcppArmadillo.h>
#include "Losses.h"

namespace rB1 {

    double ff(const colvec& L, const colvec& U, const colvec& e) {
        size_t i;
        double loss = 0.0;
        for (i = 0; i < e.size(); i++) {
            if (L[i] > e[i]) {
                loss += fabs(e[i] - L[i]) / fabs(L[i]);
                continue;
            }
            if (U[i] < e[i]) {
                loss += fabs(U[i] - e[i]) / U[i];
            }
        }
        return loss;
    }

    template <typename T> colvec ffGrd(const T& A, const colvec& e, const colvec& L, const colvec& U) {
        size_t i;
        colvec Df = zeros<colvec>(e.size());
        for (i = 0; i < e.size(); i++) {
            if (L[i] > e[i]) {
                Df[i] = 1.0 / fabs(L[i]); // d_sign(e[i] - L[i])
                continue;
            }
            if (U[i] < e[i]) {
                Df[i] = -1.0 / U[i]; // d_sign(e[i] - U[i])
            }
        }
        colvec grd = A.t() * conv_to<colvec>::from(Df);
        return grd;
    }

    template <typename T> int updategrd(const T& A, const mat& B, const colvec& s, const colvec& ee, colvec& grad, umat& ord, int j) {
        size_t i;
        bool ch = false;
        colvec Df = conv_to<colvec>::from(ee > B.col(1)) / B.col(1) - conv_to<colvec>::from(ee < B.col(0)) / abs(B.col(0));
        colvec u = conv_to<colvec>::from(s > 0) / B.col(1) - conv_to<colvec>::from(s < 0) / abs(B.col(0));
        u = Df - u;
        for (i = 0; i < u.size(); i++) {
            if (u[i] != 0.0) {
                grad -= A.row(i).t() * u[i];
                ch = true;
            }
        }
        if (ch) {
            ord = stable_sort_index(abs(grad), "descend");
            j = -1;
        }
//colvec gradstd = ffGrd(A, ee, B.col(0), B.col(1));
//cout << (gradstd - grad) - (A.t() * d) << endl;
        return j;
    }

}

#endif // RB1_HPP

