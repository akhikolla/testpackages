// #ifndef DDC_H
// #define DDC_H
// #endif
// 
// #ifndef ARMA_DONT_PRINT_ERRORS
// #define ARMA_DONT_PRINT_ERRORS
// #endif
// 
// #ifndef  ARMA_USE_CXX11
// #define ARMA_USE_CXX11
// #endif
// 
// 
// #include "RcppArmadillo.h"
// 
// 
// struct invMat {
//   arma::mat Inv;
//   arma::mat InvSqrt;
// };
// 
// invMat mpinv(arma::mat &S) {
//   invMat output;
//   arma::mat U;
//   arma::vec s;
//   arma::mat V;
//   
//   svd(U,s,V,S);
//   arma::vec sinv = 1 / s:
//   arma::vec sinvsqrt = 1 / arma::sqrt(s):
//   sinv(arma::find(arma::abs(s) < 1e-12)) = 0;
//   sinvsqrt(arma::find(arma::abs(s) < 1e-12)) = 0;
//   output.Inv =  U*diagmat(sinv)*V.t();
//   
//   return(output);
// }
// 
// arma::vec huberweights(arma::vec x, double b} {
//   arma::vec result(x.size(), arma::fill::ones);
//   outind = arma::find(arma::abs(z-mu) > 1.5);
//   result(outind) = 1.5 / arma::abs(result(outind));
//   return(result);
// }
