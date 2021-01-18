#include <RcppArmadillo.h>

//' Intraclass correlation used to assess variability of lower-order
//' relationships across higher-order processes/units.
//'
//' A function and vignettes for computing the intraclass correlation described
//' in Aguinis & Culpepper (2015). iccbeta quantifies the share of variance
//' in an outcome variable that is attributed to heterogeneity in slopes due to
//' higher-order processes/units.
//' 
//' @param X    The design `matrix` of fixed effects from a lmer model.
//' @param l2id A `vector` that identifies group membership. The vector
//'             must be coded as a sequence of integers from 1 to J,
//'             the number of groups.
//' @param T    A `matrix` of the estimated variance-covariance matrix of
//'             a lmer model fit.
//' @param vy   The variance of the outcome variable.
//' 
//' @return 
//' A `list` with:
//' 
//' - `J`
//' - `means`
//' - `XcpXc`
//' - `Nj`
//' - `rho_beta`
//' 
//' @author 
//' Steven Andrew Culpepper
//' 
//' @seealso 
//' 
//' [lme4::lmer()], [model.matrix()],
//' [lme4::VarCorr()], [RLRsim::LRTSim()],
//' [iccbeta::Hofmann], and [iccbeta::simICCdata]
//'          
//' @references
//' Aguinis, H., & Culpepper, S.A. (2015). An expanded decision making
//' procedure for examining cross-level interaction effects with multilevel
//' modeling. _Organizational Research Methods_. Available at:
//' <http://hermanaguinis.com/pubs.html>
//' 
//' @noRd
// [[Rcpp::export]]
Rcpp::List icc_beta_cpp(const arma::mat &X, const arma::vec &l2id,
                        const arma::mat &T, double vy)
{
  
  if(!X.is_finite()) {
    Rcpp::stop("`X` must not have any missing values (NA).");
  }
  
  unsigned int N = l2id.n_elem;
  unsigned int p = X.n_cols;
  unsigned int J = max(l2id);
  arma::mat SUMS = arma::zeros<arma::mat>(J, p);
  arma::mat means = arma::zeros<arma::mat>(J, p);
  arma::vec Nj = arma::zeros<arma::vec>(J);
  arma::vec ONE = arma::ones<arma::vec>(N);
  
  for (unsigned int i = 0; i < N; ++i) {
    // assume l2id goes from 0 to J-1
    Nj(l2id(i) - 1) += ONE(i);
    SUMS.row(l2id(i) - 1) += X.row(i);
  }
  
  for (unsigned int j = 0; j < J; ++j) {
    if (Nj(j) > 0.0) {
      means.row(j) = SUMS.row(j) / Nj(j);
    }
  }
  
  arma::mat XcpXc = arma::zeros<arma::mat>(p, p);
  for (unsigned int i = 0; i < N; ++i) {
    XcpXc += (X.row(i) - means.row(l2id(i) - 1)).t() *
      (X.row(i) - means.row(l2id(i) - 1));
  }
  
  double rho_beta = trace(T * XcpXc / (sum(Nj) - 1.0)) / vy;
  
  return Rcpp::List::create(Rcpp::Named("J", J),
                            Rcpp::Named("means", means),
                            Rcpp::Named("XcpXc", XcpXc),
                            Rcpp::Named("Nj", Nj),
                            Rcpp::Named("rho_beta", rho_beta));
}
