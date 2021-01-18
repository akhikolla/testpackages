#ifndef RB2_HPP
#define RB2_HPP

#include <RcppArmadillo.h>

namespace rB2 {

    double ff(const colvec& L, const colvec& U, const colvec& e) {
        size_t i;
        double tmpl, loss = 0.0;
        for (i = 0; i < e.size(); i++) {
            if (L[i] > e[i]) {
                tmpl = e[i] - L[i];
                loss += (tmpl * tmpl) / fabs(L[i]);
            }
            if (U[i] < e[i]) {
                tmpl = U[i] - e[i];
                loss += (tmpl * tmpl) / U[i];
            }
        }
        return loss;
    }

    template <typename T> colvec ffGrd(const T& A, const colvec& e, const colvec& L, const colvec& U) {
        size_t i;
        colvec Df = zeros<colvec>(e.size());
        for (i = 0; i < e.size(); i++) {
            if (L[i] > e[i]) {
                Df[i] = 2.0 * (e[i] - L[i]) / fabs(L[i]);
            }
            if (U[i] < e[i]) {
                Df[i] = 2.0 * (e[i] - U[i]) / U[i];
            }
        }
        colvec grd = - A.t() * Df;
        return grd;
    }

}

#endif // RB2_HPP
