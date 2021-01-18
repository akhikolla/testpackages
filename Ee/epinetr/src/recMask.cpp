#include <Rcpp.h>
using namespace Rcpp;

// Converts a logical vector with recombination points into a recombination mask.

// [[Rcpp::export]]
LogicalVector recMask(LogicalVector bvec)
{
  LogicalVector bvec2(bvec.size());

  bool flag = true;

  for (int i=0; i < bvec.size(); i++) {
    if (bvec[i])
      flag = !flag;

    bvec2[i] = flag;
  }

  return (bvec2);
}


// NumericVector stdnorm(int n, double sdev) {
//   NumericVector vec(n);
//
//   vec = rnorm(n);
//   vec = (vec - mean(vec)) / sd(vec) * sdev;
//
//   return(vec);
// }
//
//
// //' @export
// // [[Rcpp::export]]
// NumericVector foo(NumericVector vec, double sdev, double tolerance, int iterMax)
// {
//   int n = vec.size();
//   NumericVector best;
//
//   best = stdnorm(n, sdev);
//
//   for (int i=0; i<iterMax; i++) {
//
//   }
//
//   return(best);
// }
