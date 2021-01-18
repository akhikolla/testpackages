#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// sigma2marginal returns the marginal posterior value with respect to
// a linear regression model evaluated at sigma2.

double sigma2marginaluni(double sigma2, const arma::mat& XtX,
                      const arma::mat& SigmaBetaInv, const arma::mat& Xstar,
                      const arma::vec& Xty, const arma::vec& mu0,
                      const arma::vec& ystar) {

  double l1, l2, l3;

  int p      = XtX.n_rows;
  int nplusp = Xstar.n_rows;
  int n      = nplusp-p;

  arma::mat Vbeta_inv(p,p);

  arma::vec SigmaStar = arma::zeros<arma::vec>(n);
  arma::vec SigmaBetaInvStar = diagvec(SigmaBetaInv);

  // Create a vectorized covariance matrix
  for (int i=0; i<n; i++){
    SigmaStar(i) = 1/sigma2;
  }

  SigmaStar = join_cols(SigmaStar, SigmaBetaInvStar);

  // Compute the mean vector and covariance matrix of the regression coefs
  Vbeta_inv = XtX/sigma2 + SigmaBetaInv;
  arma::vec beta_hat  = inv(Vbeta_inv) * (Xty/sigma2 + SigmaBetaInv*mu0);

  // Compute the residual term for the quadratic
  arma::vec eta = ystar - Xstar * beta_hat;

  // Compute the log-likelihood
  l1 = -0.5 * (n+2) * log(sigma2) ;
  l2 = -0.5 * log(det(Vbeta_inv));
  l3 = -0.5 * sum( eta%eta%SigmaStar );

  // Exponentiate the log-likelihood and add a small value to avoid underflow
  double logL = l1+l2+l3; //exp(l1+l2+l3) + 0.000000001;
  return logL;
}



// [[Rcpp::export]]
SEXP sigma2marginal(int n, const arma::vec& grid,
                            const arma::mat& XtX,
                            const arma::mat& SigmaBetaInv,
                            const arma::mat& Xstar,
                            const arma::vec& Xty, const arma::vec& mu0,
                            const arma::vec& ystar) {


  arma::vec sigma2Marginal = arma::zeros<arma::vec>(n);

  for (int i=0; i<n; i++){
    sigma2Marginal(i) = sigma2marginaluni(grid(i), XtX, SigmaBetaInv,
                                       Xstar, Xty, mu0, ystar);
  }

  return Rcpp::List::create(Rcpp::Named("sigma2") = grid,
                            Rcpp::Named("logL")   = sigma2Marginal);
}


// [[Rcpp::export]]
SEXP sigma2marginalmc(int n, const arma::vec& grid,
                            const arma::mat& XtX,
                            const arma::mat& SigmaBetaInv,
                            const arma::mat& Xstar,
                            const arma::vec& Xty, const arma::vec& mu0,
                            const arma::vec& ystar) {
  
                            
  arma::vec sigma2Marginal = arma::zeros<arma::vec>(n);

  // Set-up indicator for randomly selecting rows of SigmaBetaInv and
  // thus randomly selecting values of alpha
  int m = SigmaBetaInv.n_rows;
  arma::vec Sigmarandom     = arma::randi<arma::vec>(m, distr_param(0,m-1));
  
  for (int i=0; i<n; i++){
    arma::mat SigmaBetaInvMat = arma::diagmat(SigmaBetaInv.row(Sigmarandom(i)));
    sigma2Marginal(i) = sigma2marginaluni(grid(i), XtX, SigmaBetaInvMat,
                                          Xstar, Xty, mu0, ystar);
  }

  return Rcpp::List::create(Rcpp::Named("sigma2")         = grid,
                            Rcpp::Named("logL")           = sigma2Marginal,
                            Rcpp::Named("SigmaBetaInvID") = Sigmarandom);
}





arma::mat Rmvn(int n, const arma::vec& mu_,
               const arma::mat& Sigma_) {

  int p      = Sigma_.n_rows;

  // Set-up matrix to fill with samples
  arma::mat beta(p,n);

  //Choleskey decomp of Sigma_
  arma::mat L = trans(chol(Sigma_));

  // Draw n multivariate standard normal samples
  for (int i=0; i<n; i++){
    beta.col(i) = mu_ + L*randn(p);
  }

  return beta;
}


arma::mat betaRegSamplerUni(double sigma2, const arma::mat& XtX,
                    const arma::mat& SigmaBetaInv,
                    const arma::vec& mu0, const arma::vec& Xty,
                    int nsamples) {

  //Compute inverse covariance
  arma::mat Vbeta = inv(XtX/sigma2 + SigmaBetaInv);

  //Compute mean estimate of beta
  arma::vec betamean = SigmaBetaInv*mu0 + Xty/sigma2;
            betamean = Vbeta * betamean;

  //Simulate beta from multivariate normal
  arma::mat beta = trans(Rmvn(nsamples,betamean, Vbeta));

  //Append sigma2 to the beta matrix
  arma::vec sigma2vec = arma::zeros<arma::vec>(nsamples);
            sigma2vec = sigma2vec + sigma2;

  beta = join_rows(beta, sigma2vec);

  return beta;
}


// [[Rcpp::export]]
arma::mat betaRegSampler(const arma::vec& sigma2, const arma::mat& XtX,
                    const arma::mat& SigmaBetaInv,
                    const arma::vec& mu0, const arma::vec& Xty,
                    int nsamples) {

  int p = sigma2.size();

  arma::mat Beta;
  arma::mat Beta0;

  for (int i=0; i<p; i++){
    arma::mat Beta0 = betaRegSamplerUni(sigma2(i), XtX, SigmaBetaInv, mu0,
                                        Xty, nsamples);
    Beta = join_cols(Beta, Beta0);
  }

  return Beta;
}


// [[Rcpp::export]]
arma::mat betaRegSamplermc(const arma::vec& sigma2, const arma::mat& XtX,
                    const arma::mat& SigmaBetaInv, const arma::vec& SigmaBetaInvID,
                    const arma::vec& mu0, const arma::vec& Xty,
                    int nsamples) {

  int p = sigma2.size();

  arma::mat Beta;
  arma::mat Beta0;

  for (int i=0; i<p; i++){
    arma::mat SigmaBetaInvMat = arma::diagmat(SigmaBetaInv.row(SigmaBetaInvID(i)));
    
    arma::mat Beta0 = betaRegSamplerUni(sigma2(i), XtX, SigmaBetaInvMat, mu0,
                                        Xty, nsamples);
    Beta = join_cols(Beta, Beta0);
  }

  return Beta;
}





