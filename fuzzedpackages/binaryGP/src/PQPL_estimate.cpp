#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List PQPL_estimate(arma::mat M, arma::mat X, arma::vec w, arma::vec etaTilde,
                        arma::mat R, double sigma2, int T) {

  int n = X.n_rows;
  int N = M.n_rows;
  int d = M.n_cols;
  arma::mat SigmaVinv(N, n);        // SigmaV^(-1)
  arma::mat VMinv(N, d);            // V^(-1)M
  SigmaVinv.fill(0);
  VMinv.fill(0);

  for (int t = 0 ; t < T ; t++) {
    arma::mat WBinv(n, n);
    WBinv.fill(0);
    for(int i = 0; i < n ; i++) WBinv(i,i) = 1/w(t*n+i);
    arma::mat VBinv =  inv(sigma2 * R  + WBinv);
    arma::mat SigmaVBinv = sigma2 * R * VBinv;
    arma::mat MB(n, d);
    for(int i = 0; i < n ; i++) MB.row(i) = M.row( t*n+i );
    arma::mat VMBinv = VBinv * MB;

    for(int i = 0; i < n ; i++) SigmaVinv.row( t*n+i ) = SigmaVBinv.row(i);
    for(int i = 0; i < n ; i++) VMinv.row( t*n+i ) = VMBinv.row(i);
  }

  arma::mat beta = inv(M.t() * VMinv) * VMinv.t() * etaTilde; // beta_hat
  arma::vec residual = etaTilde - M * beta;
  arma::vec zHat(N);                                          // Z_hat
  for (int t = 0 ; t < T ; t++) {
    arma::mat SigmaVBinv(n,n);
    arma::vec residualB(n);
    for(int i = 0; i < n ; i++) SigmaVBinv.row(i) = SigmaVinv.row( t*n+i );
    for(int i = 0; i < n ; i++) residualB(i) = residual( t*n+i );
    arma::vec zHatB = SigmaVBinv * residualB;
    for(int i = 0; i < n ; i++) zHat( t*n+i ) = zHatB(i);
  }

  arma::mat eta = M * beta + zHat;
  return (List::create(Named("eta") = eta, Named("beta") = beta));
}
