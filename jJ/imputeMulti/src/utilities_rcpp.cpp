#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <string>

using namespace Rcpp;
using namespace std;

//sup of L1 distance between x and y
// [[Rcpp::export]]
double supDistC (const NumericVector& x, const NumericVector& y) {
  int nx = x.size();

  double sup = -1.0;
  for (int i = 0; i < nx; i++) {
    if (abs(x[i] - y[i]) > sup) sup = abs(x[i] - y[i]);
  }
  return sup;
}

// @description Compare an array with missing values \code{marg} and an array
// with complete values \code{complete}. Return matching indices. Can compare
// either marginal-to-complete or complete-to-marginal.
// @param marg A two dimensional array with missing values
// @param complete A two dimensional array without missing values
// @param marg_to_comp Logical. Do you wish to compare marginal values to
// complete values/matches? Defaults to \code{FALSE} ie- complete values compared
// to marginal matches.
// @return A \code{list} of matches.
// [[Rcpp::export]]
List marg_comp_compare (IntegerMatrix& marg, IntegerMatrix& complete, const bool marg_to_complete= false) {

  int ncol_c = complete.ncol();
  int nrow_m = marg.nrow(), nrow_c = complete.nrow();

  if (marg_to_complete == false) { // compare rows of complete to rows of marg
    vector<vector<int> > out(nrow_c); // output vector of vectors, will be coerced to list

    for (int ci = 0; ci < nrow_c; ci++) {
      IntegerVector comp = complete.row(ci);
      for (int mj = 0; mj < nrow_m; mj++) {
        IntegerVector comp2 = marg.row(mj);
        LogicalVector eq;
        for (int k = 0; k < ncol_c; k++) {
          // check for missing values
          if (comp[k] != NA_INTEGER && comp2[k] != NA_INTEGER) {
            eq.push_back(comp[k] == comp2[k]);
          }
        }
        if(is_true(all(eq))) out[ci].push_back(mj + 1); // use R's iterator syntax (ie [1,n] vs [0, n-1]
      }
    }
    return wrap(out);
  }
  else { // compare rows of marg to rows of complete
    vector<vector<int> > out(nrow_m); // output vector of vectors, will be coerced to list

    for (int mi = 0; mi < nrow_m; mi++) {
      IntegerVector comp = marg.row(mi);
      for (int cj = 0; cj < nrow_c; cj++) {
        IntegerVector comp2 = complete.row(cj);
        LogicalVector eq;
        for (int k = 0; k < ncol_c; k++) {
          // check for missing values
          if (comp[k] != NA_INTEGER && comp2[k] != NA_INTEGER) {
            eq.push_back(comp[k] == comp2[k]);
          }
        }
        if(is_true(all(eq))) out[mi].push_back(cj + 1); // use R's iterator syntax (ie [1,n] vs [0, n-1]
      }
    }
    return wrap(out);
  }
}

