//[[Rcpp::depends(gRbase)]]
#include <RcppArmadillo.h>
#include <gRbase.h>

using namespace Rcpp;
using namespace gRbase;

//[[Rcpp::export]]
NumericVector fit2way_ (const NumericVector& tab1,
			const NumericVector& tab2,
			const CharacterVector& R,
			const CharacterVector& vn) {

  if (R.length() > 0){
    NumericVector tab3 = tab_marg_(tab1, R);
    NumericVector out = tab_perm_(tab_div_(tab_mult_(tab1, tab2), tab3), vn);
    return out;
  } else {
    double s = sum(tab1);
    NumericVector out = tab_perm_(tab_mult_(tab1, tab2), vn) ;
    out = out / s;
    return out;
  }
}


// //[[Rcpp::export]]
// NumericVector fit2way1_ (const NumericVector& tab1,
// 			 const NumericVector& tab2,
// 			 const NumericVector& tab3,
// 			 const CharacterVector& vn) {
//   //Rcout << "hej heh " << std::endl;
//   NumericVector out = tab_perm_(tab_div_(tab_mult_(tab1, tab2), tab3), vn);
//   return out;
// }

// //[[Rcpp::export]]
// NumericVector fit2way2_ (const NumericVector& tab1,
// 			 const NumericVector& tab2,
// 			 const CharacterVector& vn) {

//   double s = sum(tab1);
//   NumericVector out = tab_perm_(tab_mult_(tab1, tab2), vn) ;
//   out = out / s;
//   return out;
// }
