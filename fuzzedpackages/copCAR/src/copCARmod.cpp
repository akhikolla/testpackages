#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
using namespace std;
using namespace arma;

//src = '
//    using namespace arma;
//    mat B = Rcpp::as<mat>(B_);
//    mat M = Rcpp::as<mat>(M_);
//    colvec lam = Rcpp::as<colvec>(lam_);
//    int n = lam.n_elem;
//    colvec temp(n);
//    int k = Rcpp::as<int>(k_);
//    for (int j = 0; j < k; j++)
//    {
//    	for (int i = 0; i < n; i++)
//    	    temp[i] = pow(lam[i], j + 1);
//        M.col(j + 1) = B * temp;
//    }
//    return Rcpp::wrap(M);
//'

// [[Rcpp::export]]
arma::mat buildM_(const mat& B, int k, const colvec& eigval)
{
    int n = eigval.n_elem;
    mat M = ones<mat>(n, k + 1);
    colvec temp(n);
    for (int j = 0; j < k; j++)
    {
        for (int i = 0; i < n; i++)
            temp[i] = pow(eigval[i], j + 1);
        M.col(j + 1) = B * temp;
    }
    return M;
}

RCPP_MODULE(copCARmod)
{
    Rcpp::function("buildM_", &buildM_);
}
