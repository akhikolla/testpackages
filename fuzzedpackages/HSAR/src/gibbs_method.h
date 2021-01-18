#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma;

double HSAR_draw_rho(const mat& ,const mat& ,const mat& , const mat& ,const mat& ,const mat& , const mat& , double );

double HSAR_draw_lambda(const mat& ,const mat& ,const mat& , const mat& , double );
 
