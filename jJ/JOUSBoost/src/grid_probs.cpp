#include <Rcpp.h>
using namespace Rcpp;
#define TOL 1e-8

//' Function to compute predicted quantiles
//'
//' Find predicted quantiles given classification results at different quantiles.
//'
//' @param X Matrix of class predictions, where each column gives the predictions
//'        for a given quantile in q.
//' @param q The quantiles for which the columns of X are predictions.
//' @param delta The number of quantiles used.
//' @param median_loc Location of median quantile (0-based indexing).
//' @export
// [[Rcpp::export]]
NumericVector grid_probs(IntegerMatrix X, NumericVector q, double delta,
                       int median_loc) {

  // Check to make sure the median is located at median_loc
  if(fabs(q[median_loc] - 0.5) > TOL)
    stop("Check your median location!");
  int n_pred = X.nrow(), i, j;
  int n_q = X.ncol();
  NumericVector phat(n_pred);

  for(i = 0; i < n_pred; i++){
    if(X(i,median_loc) == 1){
      phat[i] = 1 - 1/(2*delta);
      for(j = median_loc + 1; j < n_q; j++)
        if(X(i,j) == -1){
          phat[i] = q[j] - 1/(2*delta);
          break;
        }
      } else{
        phat[i] = 1/(2*delta);
        for(j = median_loc - 1; j >= 0; j--)
          if(X(i,j) == 1){
            phat[i] = q[j] + 1/(2*delta);
            break;
          }
      }
  }

  return phat;

}

