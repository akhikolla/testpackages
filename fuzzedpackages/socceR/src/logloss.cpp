#include <Rcpp.h>
using namespace Rcpp ;

//' Computes the log loss score for a tournament prediction
//'
//' @description Compute the (weighted) rank probability score for a tournament.
//' @param m An R*T prediction matrix where the R rows represent the ordered ranks and each column is a team. Each column should sum to 1, and each row should sum to the number of teams that can attain a given rank. 
//' @param outcome A vector of length T containing the integers 1 to R giving the ranks that were obtained by each of the T teams
//' @param rankweights A vector of length R of rank weights or a single weight which will be reused for all ranks (defaults to 1)
//' @return The rank probability score. Zero means a perfect score.
//' @author Claus Ekstrom <ekstrom@@sund.ku.dk>
//' @examples
//'
//' m1 <- matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, .5, .5, 0, 0, .5, .5), 4)
//' m1 # Prediction where certain on the top ranks
//' logloss(m1, c(1, 2, 3, 4)) 
//'
//' @export
// [[Rcpp::export]]
double logloss(const NumericMatrix& m, NumericVector outcome, NumericVector rankweights=1) {
  double result = 0.0;

  // Expand weight if it was just a single number
  if (rankweights.size() != m.nrow())
    rankweights = NumericVector(m.nrow(), rankweights(0));

  // Should the weight be normalized?
  // This isn't implemented in the code here ... yet

  // Iterate over teams (columns)
  for (int j = 0; j < m.ncol(); ++j) {
    result += rankweights(outcome(j)-1)*log(m(outcome(j)-1, j));    
  }
  return -result;
}
