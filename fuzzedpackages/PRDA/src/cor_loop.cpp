// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
List cor_loop(int n, String alternative, int B, arma::mat Eigen_matrix) {

  NumericVector estimate(B);
  NumericVector pval(B);
  double tstat = 0;
  int n1 = n - 1;
  double df = n - 2;

  for(int i = 0; i < B; i++){
    // sample observations
    NumericVector v = rnorm(2 * n);
    v.attr("dim") = Dimension(n, 2);

    arma::mat out = as<arma::mat>(v);
    out = Eigen_matrix * out.t();

    arma::vec x = out.rows(0, 0).t();
    arma::vec y = out.rows(1, 1).t();

    // compute cor
    double mx = mean(x);
    double my = mean(y);

    double cor = 0;
    double sdx = arma::stddev(x);
    double sdy = arma::stddev(y);
    for (int j = 0; j < n; j++){
      cor += (x[j] - mx) * (y[j] - my);
    }

    cor /= n1;
    cor /= (sdx * sdy);

    // compute test
    tstat = sqrt(df) * cor / sqrt(1 - pow(cor, 2));

    if (alternative == "two_sided"){
      pval[i] = 2*(Rf_pt(std::abs(tstat), df, false, false));
    } else if (alternative == "greater") {
      pval[i] = Rf_pt(tstat, df, false, false);
    } else {
      pval[i] = 1-Rf_pt(tstat, df, false, false);
    }

    estimate[i] = cor;

  }

  return List::create(Named("p_value")=pval,
                      Named("estimate")=estimate);
}


//


