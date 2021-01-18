// [[Rcpp::depends(RcppArmadillo)]]
#pragma GCC optimize ("O3")
#include <RcppArmadillo.h>
#include <Rcpp.h>

/*
	Auxiliary function for the Predictive Mean Matching in BBPMMrow
 Version:             0.2
 Date:         2014-12-17
	Author:             T.S.
 Further infos, references and credits:
  See for Rcpp: Eddelbuettel, D. and Francois, R. (2011) Rcpp: Seamless R and C++ Integration.
                Journal of Statistical Software, Vol. 40, No. 8, pp. 1--18. URL http://www.jstatsoft.org/v40/i08/.

  See for RcppArmadillo: Eddelbuettel, D. and Sanderson, C. (2014) RcppArmadillo: Accelerating R with high-performance C++ linear algebra.
                         Computational Statistics and Data Analysis, Vol. 71, March 2014, pp. 1054--1063.
 License:  (GPL >= 2)
*/

using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec PMMC(arma::colvec vc_yhat, arma::mat mtx_weight, arma::mat mtx_yobs){
	try{
    //Def.& Dec.
    arma::rowvec vr_ans(mtx_yobs.n_cols);
	
	for(unsigned int i = 0; i < mtx_yobs.n_cols; i++) {
		vr_ans(i) = arma::as_scalar(arma::trans(vc_yhat - mtx_yobs.col(i))*mtx_weight*(vc_yhat - mtx_yobs.col(i)));
    }
              
    //Return result
    return vr_ans;
} catch (...){
   ::Rf_error("Unknown C++ Exception Error.");
}
}

RcppExport SEXP PMMC(SEXP vc_yhatSEXP, SEXP mtx_weightSEXP, SEXP mtx_yobsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::colvec >::type vc_yhat(vc_yhatSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type mtx_weight(mtx_weightSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type mtx_yobs(mtx_yobsSEXP );
        arma::rowvec __result = PMMC(vc_yhat, mtx_weight, mtx_yobs);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


