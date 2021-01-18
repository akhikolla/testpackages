// [[Rcpp::depends(RcppArmadillo)]]
#include "post.h"




//' Helper function to sample covariance
//' 
//' This function is called by \code{sur_sample_cov_cpp}.
//' It samples the covariance matrix of a SUR
//'
//' @param Y A \code{matrix}, each column a \code{vector} of responses
//' @param Xlist A \code{list}, each element a design \code{matrix}
//' @param n Integer giving number of observations
//' @param J Integer giving number of endpoints
//' @param pj A \code{vector} giving number of covariates per endpoint
//' @param sigma11 A scalar giving a draw for the (1,1) component of the covariance matrix
//' @param r1 A \code{vector} of residuals for the first endpoint's regression
// [[Rcpp::export]]
arma::mat sur_sample_cov_helper_cpp(
    arma::mat const& Y,
    List const& Xlist, 
    int const& n, 
    int const& J, 
    arma::vec const& pj,
    double const& sigma11,
    arma::vec const& r1
) {
  arma::mat Rj = r1;                   // matrix of residuals, initialized to be first residual given as param
  arma::mat Sigma(J, J, fill::zeros);   // covariance matrix to sample
  Sigma(0,0) = sigma11;         // initialize Sigma[1,1] element
  
  for ( int j = 1; j < J; j++ ) {
    
    // Get yj, Xj, and augmented design matrix
    //   Zj = (r1, ..., r_{j-1}, Xj)
    arma::vec yj = Y.col(j);
    arma::mat Xj = Xlist[j];
    arma::mat Zj = join_rows(Rj, Xj);
    
    // Get MLE covariance and mean estimates for yj ~ Zj
    arma::mat covbhatj = inv_sympd(Zj.t() * Zj);                         // (Zj' Zj)^{-1}
    arma::vec bhatj = covbhatj * (Zj.t() * yj);                        // (Zj' Zj)^{-1} Zj' y
    
    // Sample omegasqj ~ InvGamma( 0.5 * nuj, 0.5 * SSEj ) See (14) in Zellner.
    int nuj = n - Zj.n_cols;
    double ssej = pow( norm(yj - Zj * bhatj), 2);              // squared norm of residual vector
    double omegasqj = 1.0 / R::rgamma(nuj / 2.0, 2.0 / ssej );  // uses rate paramaterization instead of scale
    
    // Sample bj = N(bhatj, omegasqj * (Zj' Zj)^(-1) ). See (13) in Zellner
    arma::vec bj = mvnrnd(bhatj, omegasqj * covbhatj);
    
    // Update residual matrix w/ elements of bj corresponding to betaj (Xj)
    arma::vec betaj = bj.tail(pj(j));              // last pj elements of bj
    Rj = join_rows(Rj, yj - Xj * betaj);
    
    // Update sigma matrix according to eq (9) of Zellner & Ando
    if( j == 1 ) {
      Sigma(1,1) = omegasqj + pow(bj(0), 2) * Sigma(0,0);   // update Sigma[2,2] = omegasq + rho1^2 * Sigma[1,1]
      Sigma(1,0) = Sigma(0,1) = bj(0) * Sigma(0,0);
      continue;
    }
    
    if ( j > 1 ) {
      arma::vec rhoj = bj.head(j);    // first j elements of bj
      
      // Formula for Sigma[j,i] is sum_{k=1}^{j-1} \rho_{jk} \sigma_{ki}
      //   
      //   This is just the dot product t(Sigma[, i] ) %*% rhoj, i = 1, 2, ..., j-1
      //   We can get all upper triangular elements simultaneously by performing 
      //   Sigmasubj %*% rhoj where Sigmasubj is (j-1) x (j-1) submatrix of Sigma
      //   
      arma::mat Sigmasubj = Sigma.submat( 0, 0, j-1, j-1 ); // (firstrow, firstcol, lastrow, lastcol)
      arma::vec uppertrij = Sigmasubj * rhoj;
      Sigma.submat( 0, j, j-1, j ) = uppertrij;   // (firstrow, firstcol, lastrow, lastcol)
      
      //   Formula for Sigma[j,j] is
      //   Sigma[j,j] = omega_j^2 
      //                  + sum_{k=1}^{j-1} \rho_{jk}^2 \sigma_{kk} 
      //                  + \sum_{k \ne l}^{j-1} \rho_{jk} \rho_{jl} \sigma_{lk}
      //  Upon closer inspection, the two sums combine to form the qudratic form
      //               \rho_{j}^T Sigma_{j-1} \rho_{j}
      //  where Sigma_{j-1} is the matrix from the previous iteration
      Sigma(j,j) = omegasqj + dot(rhoj, Sigmasubj * rhoj);   // rhoj' * Sigmasubj * rhoj
      
      // fill in the lower-triangular elements
      Sigma = symmatu(Sigma);                                
      
    }  // end if j > 1 statement
  }    // end j loop
  return( Sigma );
}




//' Sample from SUR via Direct Monte Carlo (C++ version)
//' 
//' C++ implementation of Zellner and Ando (2010) Direct Monte Carlo
//' method for sampling from the posterior of a Bayesian SUR
//'
//' @param Y \code{matrix} \eqn{(y_1, \ldots y_J)}
//' @param Xlist A \code{list}, each element a design \code{matrix}
//' @param y \code{vector} of responses
//' @param X design \code{matrix}
//' @param XtX \code{matrix} giving \code{crossprod(cbind(X1, ..., XJ))}
//' @param pj \code{vector} giving number of covariates per endpoint
//' @param M An integer giving the number of desired samples
// [[Rcpp::export]]
List sur_sample_cpp (
    arma::mat const& Y,
    List const& Xlist,
    arma::vec const& y,
    arma::mat const& X,
    arma::mat const& XtX,
    arma::vec const& pj,
    int const& M
) {
  // Empty list to store covariances
  List Sigmalist(M);
  
  // Empty matrix to store drawn betas
  arma::mat betasample(sum(pj), M, fill::none);
  
  // Get vars for first endpoint
  arma::mat const& X1 = Xlist[0];
  arma::vec const& y1 = Y.col(0);
  
  // Get some parameters needed for loop and calculations
  int const& J = Y.n_cols;
  int const& n = y1.size();
  
  
  // get (X1'X1)^(-1) and bhat1 = MLE for y1 ~ X1
  arma::mat const& covbhat1 = inv_sympd(X1.t() * X1);
  arma::vec const& bhat1 = covbhat1 * (X1.t() * y1);
  
  // Since Z1 = X1 (regression never changes),
  // we can sample all Sigma[1,1] values simultaneously,
  // avoiding loops
  
  // Sample Sigma11 ~ InvGamma(0.5 * nu1, 0.5 * sse1);
  // Get parameters for gamma draw--sse1 and nu1
  double const& sse1 = pow( norm(y1 - X1 * bhat1), 2);
  int const& nu1 = n - pj(0);
  
  // Must sample from gamma distribution first; then take inverse
  arma::vec const& omegasqinv = randg( M, distr_param( nu1 / 2.0, 2.0 / sse1 ) );
  arma::vec const& sigma11 = pow(omegasqinv, -1.0);
  
  // Sample b1 = beta1 ~ N(bhat1, omegasq[m] * (X1'X1)^(-1) )
  // Get M x p1 matrix of N(0,1) samples
  arma::mat const& Q = randn(M, pj(0));
  
  // repeat t(bhat) matrix M times to get matrix of means;
  // change replicate sigma11 vector into M x p0 matrix and 
  // square root. Then element-wise multiplication of Q * chol(covbhat1)
  // where Q * chol(covbhat1) gives M x p1 matrix of N(0, (X1'X1)^(-1) )
  arma::mat const& b1 = repmat(bhat1, 1, M).t()
    +  repmat( pow(sigma11, 0.5), 1, pj(0) ) %  Q * chol(covbhat1) ;
  
  // Get n x M matrix of sampled residuals for first endpoint
  arma::mat const& r1 = repmat( y1, 1, M) - X1 * b1.t();
  
  // Now, loop through M = numsamples
  for ( int m = 0; m < M; m++ ) {
    arma::vec r1_m = r1.col(m);                      // get m-th residual vector
    double sigma11_m = sigma11(m);             // get m-th Sigma[1,1] sample
    
    // call helper function to get sampled covariance matrix; store in the list
    arma::mat Sigma_m = sur_sample_cov_helper_cpp( Y, Xlist, n, J, pj, sigma11_m, r1_m );
    Sigmalist[m] = Sigma_m;
    
    // Sample beta | Sigma
    
    // Get inverse of sampled matrix
    arma::mat Sigmainv = inv_sympd(Sigma_m);
    
    // Get [ t(X) %*% kron(Sigma, I) %*% X ]^{-1} matrix 
    arma::mat betacov = inv_sympd( fastKronEye_crossprod(XtX, Sigmainv, pj, n, J) );
    
    // Get t(X) %*% kron(Sigma, I) %*% y matrix
    arma::vec betamean = betacov * ( X.t() * fastKronEye_Y(Sigmainv, Y, n, J) );
    
    // beta ~ N(betamean, betacov)
    arma::vec beta_m = mvnrnd(betamean, betacov);
    betasample.col(m) = beta_m;
  }
  
  List res = List::create(
    Named("betadraw") = betasample.t(),
    Named("covlist") = Sigmalist
  );
  
  return( res );
}

