#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List clm(SEXP Xarg, SEXP yarg) {

    // déclaration des variables
    Rcpp::NumericMatrix X1(Xarg);
    Rcpp::NumericVector y1(yarg);
    int n = X1.nrow(), k = X1.ncol();
    arma::mat X(X1.begin(), n, k, false);
    arma::colvec y(y1.begin(), y1.size(), false);
    int df = n - k;

    // résolution du modèle y ~ X
    arma::colvec coef = arma::solve(X, y);

    // calcul des estimations
    arma::colvec fitted = X*coef;

    // calcul des residus
    arma::colvec res  = y - fitted;

    // écart-type des coefficients
    double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/df;
    arma::colvec sderr = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));

    // résultats sous forme de liste nommée
    return Rcpp::List::create(Rcpp::Named("coefficients") =coef,
                              Rcpp::Named("stderr")       =sderr,
                              Rcpp::Named("df")           =df,
                              Rcpp::Named("fitted.values")=fitted,
                              Rcpp::Named("res")          =res);
}
