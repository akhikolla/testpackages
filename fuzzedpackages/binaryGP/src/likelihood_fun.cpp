#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double likelihood_fun(arma::mat M, arma::mat X, arma::vec w, arma::vec etaTilde,
                   arma::mat R, double sigma2, int T) {

  int n = X.n_rows;
  int N = M.n_rows;
  int d = M.n_cols;
  arma::mat eGammae1(1,1);          // eta'V^(-1)eta
  eGammae1.fill(0);
  double logDetV = 0;               // log(det(V)), which can be computationally simplied by matrix determinant lemma
  arma::mat VMinv(N, d);            // V^(-1)M
  VMinv.fill(0);

  for (int t = 0 ; t < T ; t++) {
    arma::mat WBinv(n, n);
    WBinv.fill(0);
    for(int i = 0; i < n ; i++) WBinv(i,i) = 1/w(t*n+i);
    arma::mat WB(n, n);
    WB.fill(0);
    for(int i = 0; i < n ; i++) WB(i,i) = w(t*n+i);
    arma::mat VBinv =  inv(sigma2 * R  + WBinv);
    arma::mat MB(n, d);
    for(int i = 0; i < n ; i++) MB.row(i) = M.row( t*n+i );
    arma::mat VMBinv = VBinv * MB;

    for(int i = 0; i < n ; i++) VMinv.row( t*n+i ) = VMBinv.row(i);
    arma::vec etaTildeB(n);
    for(int i = 0; i < n ; i++) etaTildeB(i) = etaTilde( t*n+i );
    eGammae1 = eGammae1 + etaTildeB.t() * VBinv * etaTildeB;
    logDetV = logDetV + log(det(arma::eye(n,n) + WB * (sigma2 * R))) - sum(log(diagvec(WB)));
  }

  arma::mat MtVinvM = M.t() * VMinv;
  arma::mat eGammae(1,1);
  eGammae = eGammae1 - etaTilde.t() * VMinv * inv(MtVinvM) * VMinv.t() * etaTilde;

  return(logDetV + log(det(MtVinvM)) + eGammae(0,0));
}
