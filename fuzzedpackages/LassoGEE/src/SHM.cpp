#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::List SHM(arma::mat X, arma::vec y, arma::vec mu,
               arma::vec mu_eta, arma::vec vari, arma::vec nt,
               arma::vec index, arma::cube Rhat, int N, double fihat) {

  int nx = X.n_cols;
  arma::mat eachS(nx, N, arma::fill::zeros);
  arma::cube eachH(nx, nx, N, arma::fill::zeros);
  arma::cube eachM(nx, nx, N, arma::fill::zeros);

  arma::mat S(eachS.n_rows, 1, arma::fill::zeros);
  arma::mat H(eachH.n_rows, eachH.n_cols, arma::fill::zeros);
  arma::mat M(eachM.n_rows, eachM.n_cols, arma::fill::zeros);

  for (int i = 0; i < N; ++i) {
    arma::mat ym(nt[i], 1, arma::fill::zeros);
    arma::mat bigD(nt[i], nx, arma::fill::zeros);
    arma::mat bigA(nt[i], nt[i], arma::fill::zeros);

    for (int j = 0; j < nt[i]; ++j) {
      ym[j] = y[j + index[i]] - mu[j + index[i]];
      bigA(j, j) = vari[j + index[i]];
      for (int k = 0; k < nx; ++k) {
        bigD(j, k) = mu_eta[j + index[i]] * X(j + index[i], k);
      }
    }

    arma::mat Rtmp = Rhat.slice(i);
    arma::mat bigV = sqrt(bigA) * Rtmp.submat(0, 0, nt[i]-1, nt[i]-1) * sqrt(bigA);
    arma::mat sum200 = bigD.t() * inv(bigV) * ym;
    eachS.col(i) = sum200;
    arma::mat sum300 = bigD.t() * inv(bigV) * bigD;
    eachH.slice(i) = sum300;
    arma::mat sum400 = bigD.t() * inv(bigV) * ym * ym.t() * inv(bigV) * bigD;
    eachM.slice(i) = sum400;

    H += eachH.slice(i);
    M += eachM.slice(i);

  }


  for (int i = 0; i < nx; ++i) {
    S[i] = sum(eachS.row(i));
  }



  return Rcpp::List::create(Rcpp::Named("Smat") = eachS*fihat,
                            Rcpp::Named("S") = S*fihat,
                            Rcpp::Named("H") = H*fihat,
                            Rcpp::Named("M") = M*fihat);
}


