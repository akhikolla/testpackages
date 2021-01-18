// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::interfaces(r, cpp)]]

// The equivalent of the R command rowSums(x[, -j, drop = FALSE])
// that is, row sums of a matrix whose jth column has been removed

// [[Rcpp::export]]
arma::colvec arma_rowSums_minus_col(const arma::mat& x, const int & j) {
  return arma::sum(x, 1) - x.col(j - 1) ;
}

// The equivalent of the R command ifelse(x == 0, constant, log(x))

// [[Rcpp::export]]
arma::mat cpp_log0const(const arma::mat& x, const double& constant) {
  int nr = x.n_rows ;
  int nc = x.n_cols ;
  arma::mat res(nr, nc) ;
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      if (x(i,j) > 0) {
        res(i, j) = log(x(i,j)) ;
      } else {
        res(i, j) = constant ;
      }
    }
  }
  return res ;
}

// The equivalent of the R command apply(x, 2, function(x) sum(x ^ 2))
// that is, columnwise mean squares

// [[Rcpp::export]]
arma::rowvec cpp_col_ms(arma::mat const& x){
  return arma::sum(x % x, 0) / x.n_rows ;
}

// Function to estimate sigma2hat_dj

// [[Rcpp::export]]
Rcpp::List cpp_sigma2hat_dj(const Rcpp::List& all_max,
                            const int& b, const int& kn, const int& m,
                            const String& bias_adjust, const String& which_dj)
  {
  arma::vec ys = all_max["ys"] ;
  arma::vec xs = all_max["xs"] ;
  arma::mat yd = all_max["yd"] ;
  arma::mat xd = all_max["xd"] ;
  int ktop = xd.n_cols ;
  arma::mat nums_mat(kn, kn), FhatjMni(kn, kn) ;
  arma::colvec Fhaty(kn), sigma2dj_for_sl(2), sigma2dj(2), theta_dj(2) ;
  arma::rowvec That(kn), ThatN(kn), sigmahat2_dj(kn), sigmahat2_djN(kn) ;
  arma::mat Nhat(kn, ktop), Zhat(kn, ktop), ZhatN(kn, ktop) ;
  arma::mat Usum(ktop, kn), UsumN(ktop, kn), Bhat(ktop, kn), BhatN(ktop, kn) ;
  arma::mat data_dj(kn, 2) ;
  double constant =   -log(m - b + kn + 0.0) ;
  // Loop over the required columns of yd and xd
  for (int k = 0; k < ktop; k++) {
    arma::vec y = yd.col(k) ;
    arma::vec x = xd.col(k) ;
    for (int i = 0; i < kn ; i++) {
      arma::vec subx = x.subvec(i * b, (i + 1) * b - 1) ;
      for (int j = 0; j < kn; j++) {
        nums_mat(j, i) = sum(subx <= y(j)) ;
      }
    }
    Fhaty = arma::sum(nums_mat, 1) / m ;
    for (int i = 0; i < kn ; i++) {
      FhatjMni.col(i) = arma_rowSums_minus_col(nums_mat, i + 1) / (m - b);
    }
    Nhat.col(k) = Fhaty ;
    arma::rowvec pjn = arma::mean(FhatjMni, 0) ;
    Usum.row(k) = b * (1.0 - arma::mean(FhatjMni, 0)) ;
    UsumN.row(k) = -b * arma::mean(cpp_log0const(FhatjMni, constant), 0) ;
  }
  // BB2018
  Zhat = b * (1.0 - Nhat) ;
  That = arma::sum(Zhat, 0) / Zhat.n_rows;
  // N2015
  ZhatN = -b * log(Nhat) ;
  ThatN = arma::sum(ZhatN, 0) / ZhatN.n_rows;
  for (int i = 0; i < kn ; i++) {
    Usum.col(i) = kn * trans(That) - (kn - 1.0) * Usum.col(i) ;
    UsumN.col(i) = kn * trans(ThatN) - (kn - 1.0) * UsumN.col(i) ;
  }
  // The dimensions need fiddling.  Loop as above.
  for (int i = 0; i < kn ; i++) {
    Bhat.col(i) = trans(Zhat.row(i)) + Usum.col(i) - 2.0 * trans(That) ;
    BhatN.col(i) = trans(ZhatN.row(i)) + UsumN.col(i) - 2.0 * trans(ThatN) ;
  }
  for (int i = 0; i < ktop ; i++) {
    BhatN.row(i) = BhatN.row(i) - sum(BhatN.row(i)) / kn ;
  }
  sigmahat2_dj = cpp_col_ms(trans(Bhat)) ;
  sigmahat2_djN = cpp_col_ms(trans(BhatN)) ;
  sigma2dj_for_sl(0) = sum(sigmahat2_djN) / ktop ;
  sigma2dj_for_sl(1) = sum(sigmahat2_dj) / ktop ;
  if (bias_adjust == "N") {
    Nhat = (m * Nhat - b) / (m - b) ;
    That = arma::sum(b * (1.0 - Nhat), 0) / Nhat.n_rows;
    ThatN = arma::sum(-b * cpp_log0const(Nhat, constant), 0) / Nhat.n_rows;
  }
  // For disjoint maxima pick either the first or last value, based on which_dj
  if (which_dj == "first") {
    theta_dj(0) = 1.0 / ThatN(0)  ;
    theta_dj(1) = 1.0 / That(0)  ;
    sigma2dj(0) = sigmahat2_djN(0) ;
    sigma2dj(1) = sigmahat2_dj(0) ;
    data_dj.col(0) = -b * log(Nhat.col(0)) ;
    data_dj.col(1) = b * (1.0 - Nhat.col(0)) ;
  } else {
    theta_dj(0) = 1.0 / ThatN(ktop - 1)  ;
    theta_dj(1) = 1.0 / That(ktop - 1)  ;
    sigma2dj(0) = sigmahat2_djN(ktop - 1) ;
    sigma2dj(1) = sigmahat2_dj(ktop - 1) ;
    data_dj.col(0) = -b * log(Nhat.col(ktop - 1)) ;
    data_dj.col(1) = b * (1.0 - Nhat.col(ktop - 1)) ;
  }
  return List::create(Named("sigma2dj") = sigma2dj,
                      Named("sigma2dj_for_sl") = sigma2dj_for_sl,
                      Named("theta_dj") = theta_dj,
                      Named("data_dj") = data_dj) ;
}
