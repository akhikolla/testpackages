#include <Rcpp.h>
using namespace Rcpp ;

//' Computes the rank probability score for a tournament
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
//' trps(m1, c(1, 2, 3, 4)) 
//'
//' @export
// [[Rcpp::export]]
double trps(const NumericMatrix& m, NumericVector outcome, NumericVector rankweights=1) {
  double result = 0.0;
  double cumsum = 0.0;

  // Various sanity checks
  if (outcome.size() != m.ncol()) {
    stop("The number of teams (columns in m) must match the number of teams ranks (length of outcome)");    
  }

  if (max(outcome) != m.nrow()) {
    stop("The largest rank in the outcome must match the number of possible ranks (row in m)");    
  }

  // Expand weight if it was just a single number
  if (rankweights.size() != m.nrow())
    rankweights = NumericVector(m.nrow(), rankweights(0));

  if (m.nrow()==1)
    return(0.0);

  // Should the weight be normalized?
  // This isn't implemented in the code here ... yet

  int CO = 0;
  // Iterate over teams (columns)
  for (int j = 0; j < m.ncol(); ++j) {
    CO = 0;
    cumsum = m(0, j);
    if (outcome(j)==1) {
      CO = 1;
    }
    result += rankweights(0)*(cumsum - CO)*(cumsum - CO);
    //        Rcout << result << std::endl;
    // We have already taken the first row above so start at second row
    for (int i = 1; i < m.nrow()-1; ++i) {
      if (outcome(j) == (i+1)) {
	CO = 1;
      }
      cumsum += m(i, j);
      result += rankweights(i)*(cumsum - CO)*(cumsum - CO);
      //           Rcout << "   " << result << std::endl;
    }
  }
  return result/(m.rows()-1)/(m.cols());
}
