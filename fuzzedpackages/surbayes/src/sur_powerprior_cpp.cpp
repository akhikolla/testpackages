// [[Rcpp::depends(RcppArmadillo)]]
#include "post.h"



//' Sample Sigma via Gibbs for SUR model
//' 
//' This is a c++ implementation of sampling Sigma via Gibbs in SUR model--inverse Wishart
//'
//' @param nu degrees of freedom
//' @param V scale matrix
//' @param p dimension of covariance matrix
//' @return sampled covariance matrix
// [[Rcpp::export]]
arma::mat sample_sigma(double const& nu, arma::mat const& V, int const& p) {
  
  // Create empty matrix, fill with 0s
  arma::mat Q = arma::mat(p, p, fill::zeros);
  
  
  //
  // Fill upper-triangular elements of Q
  //
  
  // make empty matrices
  arma::mat X(p, p, fill::zeros);
  arma::mat Z(p, p, fill::zeros);
  
  // fill matrices with integers
  arma::vec idx = arma::linspace<mat>(1, p, p);   // = 1:p
  X.each_col() += idx;                            // cumulatively adds 1:p to each column of X
  Z.each_row() += trans(idx);                     // cumsum 1:p to each row of Z
  
  // Sample--diagonal elements of Q are sqrt of gamma
  
  // for (int i = 0; i < p; i++ ) {
  //   Q(i,i) = sqrt( rchisq(1, nu - i)[0] );
  // }
  
  arma::vec gamma_sample = arma::sqrt(
    arma::chi2rnd(
      nu - linspace(0, p-1, p)
    )
  );
  
  // off-diagonal are N(0,1). There are p * (p-1) / 2 off diag eleemnts
  arma::vec normal_sample = arma::randn( p * (p - 1) / 2 );
  
  // assign upper triangular elements
  Q.elem( find( Z == X ) ) = gamma_sample;
  Q.elem( find(Z < X) ) = normal_sample;
  
  // Cholesky decomposition of Wishart draw
  arma::mat C = Q.t() * arma::chol(V);        // upper triangular matrix
  
  // Cholesky decomposition of inverse--which is Sigma since Sigmainv ~ Wish; Sigma ~ Wish
  mat Cinv = solve( trimatu(C), eye(p, p) );   // solve ( trimatu ( X ) ) is faster than solve since upper tri mtx
  
  // return full matrix
  return Cinv * Cinv.t();
} 







//' Power Prior Gibbs sampling
//' 
//' This is a c++ implementation of Gibbs sampling SUR model with power prior
//'
//' @param Sigma initial value for covariance matrix
//' @param M number of samples
//' @param X design matrix for current data
//' @param X0 design matrix for historical data
//' @param XtX matrix that is \code{crossprod(cbind(X1, ..., XJ))}
//' @param X0tX0 matrix that is \code{crossprod(cbind(X01, ..., X0J))}
//' @param Y future response as matrix (Y1, ..., YJ)
//' @param Y0 historical response as matrix (Y01, ..., Y0J)
//' @param y future response as vector
//' @param y0 historical response as vector
//' @param a0 power prior parameter
//' @param pvec \code{vector} giving number of covariates per endpoint
//' @param burnin Burn-in parameter
//' @param thin Thin parameter
//' @return sampled covariance matrix
// [[Rcpp::export]]
List sur_sample_gibbs_cpp (
    arma::mat Sigma,
    int const& M,
    arma::mat const& X,
    arma::mat const& X0,
    arma::mat const& XtX,
    arma::mat const& X0tX0,
    arma::mat const& Y,
    arma::mat const& Y0,
    arma::vec const& y,
    arma::vec const& y0,
    double const& a0,
    arma::vec const& pvec,
    int burnin,
    int thin
) {
  int const& n = Y.n_rows;
  int const& n0 = Y0.n_rows;
  int const& J = Sigma.n_rows;
  int const& p = sum(pvec);
  
  // Create empty list, matrix, and list for results
  int lastiter = M - (M % thin);                            // last iteration of sampler
  int numstore = (int) floor( 1.0 * (M - burnin) / thin );        // final num of samples after burn and thin
  
  arma::mat betasample(p, numstore, fill::zeros);
  List sigmalist(numstore);
  
  for( int m = 0; m < lastiter; m++ ) {
    // 
    // Draw Beta | Sigma
    // 
    
    // Get Sigma^{-1}
    arma::mat Sigmainv = inv_sympd(Sigma);
    
    // Get mean and variance for posterior for new beta
    // beta_cov = (Cov for new posterior + Cov for power prior)^{-1}
    arma::mat beta_cov = arma::inv_sympd(
      fastKronEye_crossprod(XtX, Sigmainv, pvec, n, J) + 
        a0 * fastKronEye_crossprod(X0tX0, Sigmainv, pvec, n0, J )
    );
    
    arma::vec rhs1 =  X.t() * fastKronEye_Y( Sigmainv, Y, n, J );
    arma::vec rhs0 = X0.t() * fastKronEye_Y(Sigmainv, Y0, n0, J);
    
    arma::vec beta_mean = beta_cov * ( 
      X.t() * fastKronEye_Y( Sigmainv, Y, n, J ) +
        + X0.t() * fastKronEye_Y( Sigmainv, Y0, n0, J) 
    );
    
    // Draw new beta based on multivariate normal
    arma::vec beta = mvnrnd(beta_mean, beta_cov);
    
    // 
    // Draw Sigma | beta
    // 
    // Get matrix of residuals
    arma::mat resid0 = y0 - X0 * beta;      // residual vector for historical
    resid0.reshape(n0, J);                  // convert to n0 x J matrix of residuals (r_10, ..., r_J0)
    
    arma::mat resid = y - X * beta;         // residual vector for current data
    resid.reshape(n, J);                    // convert to n x J matrix of residuals (r_1, ..., r_J)
    
    arma::mat R0 = resid0.t() * resid0;     // JxJ matrix of dot product residuals for historical
    arma::mat R = resid.t() * resid;        // JxJ matrix of dot product residuals for current data
    
    arma::mat V = inv_sympd(R + a0 * R0);     // scale matrix for wishart
    double nu = n + a0 * n0;                  // degrees of freedom for wishart
    
    arma::mat Sigma = sample_sigma(nu, V, J);    // Call function that samples sigma from inv. wishart
    
    // Store after burn-in and thin
    if( m < burnin ) {
      continue;
    }
    
    int i = m - burnin;                      // iteration number post burn-in; i = 0 when m = burnin + 1
    if( ( i + 1 ) % thin == 0 ) {
      int j = (i + 1) / thin - 1;
      betasample.col(j) = beta;              // index of stored samples; j = 0 when m = burnin + 1 + thin
      sigmalist[j] =  Sigma;
    }
  }
  return List::create(
    Named("betasample") = betasample.t(),
    Named("Sigmalist") = sigmalist
  );
}


