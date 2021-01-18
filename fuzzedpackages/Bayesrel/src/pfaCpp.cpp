//
//  pfaArma.cpp
//
//  Created by Julius Pfadt on 10.08.20.
//

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
//[[Rcpp::export]]
List pfaArma(const arma::mat& X) {
    int k = X.n_cols;
    mat R = X;
    vec sds = 1/sqrt(X.diag());
    mat Xcor = diagmat(sds) * X * diagmat(sds);
    Xcor.diag().ones();
    mat XCorInv = inv_sympd(Xcor);
    vec smc = 1 - 1 / XCorInv.diag();
    R.diag() = smc;
    double h2 = sum(smc);
    double error = h2;
    int i = 1;
    vec eigval;
    mat eigvec;
    vec lambda;
    mat R_mod;
    while (error > .001 || i == 1) {
        eig_sym(eigval, eigvec, R);
        lambda = eigvec.col(k-1) * sqrt(eigval(k-1));
        R_mod = lambda * lambda.t();
        double h2_new = trace(R_mod);
        error = std::fabs(h2 - h2_new);
        h2 = h2_new;
        R.diag() = R_mod.diag();
        i += 1;
        if (i > 50) {
            error = 0;
        }
    }
    mat E = X - R_mod;
    vec ee = E.diag();
    List LL = List::create(Named("loadings")=lambda, _("err_var")=ee);
    return LL;
}
