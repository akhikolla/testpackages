#include <iostream>
#include <fstream>
#include <RcppArmadillo.h>
#include <stdlib.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends("RcppArmadillo")]]

//' Base function for computing SpaCC solution for single regularization value.
//' @param X A subject (n) by probe (p) data matrix
//' @param w A vector weights for adjacent probes. Should have length nprobes -1
//' @param gamma A scalar value for the regularization parameter
//' @param nu A scalar value for the step size in AMA algorithm
//' @param verbose Logical value whether progress should be printed
//' @param tol A scalar value for convergence tolerance.
//' @param maxiter Maximum number of iterations
//' @param Uinit A matrix used for warm starts with U
//' @param Vinit A matrix used for warm start with V
//' @param Laminit A matrix used for warm starts with Lam
//' @return An RcppArmadillo field object. Has three components, each holds the U,V, and Lam matrix for the current regularization
//' @export
// [[Rcpp::export]]
arma::field<arma::mat> SpaCC(arma::mat X, const arma::vec& w, const double& gamma, const double& nu, const bool& verbose, const double& tol, const int& maxiter, const arma::mat& Uinit, const arma::mat& Vinit, const arma::mat& Laminit){
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int wlen = w.n_elem;
  if( wlen != m -1) {
    throw std::invalid_argument( "incorrect weight length" );
  }
  if( (Uinit.n_rows != n) || (Uinit.n_cols != m) ){
    throw std::invalid_argument( "incorrect Uinit dimensions" );
  }
  if( (Vinit.n_rows != n) || (Vinit.n_cols != m-1) ){
    throw std::invalid_argument( "incorrect Vinit dimensions" );
  }
  if( (Laminit.n_rows != n) || (Laminit.n_cols != m-1) ){
    throw std::invalid_argument( "incorrect Laminit dimensions" );
  }

  double err = 1;
  double ferrnew = 1;
  double ferrold = 1;
  arma::mat tmp = arma::zeros<arma::mat>(n,1);
  arma::vec tmpmax = arma::zeros<arma::vec>(2);
  arma::vec sigma = (gamma/nu)*w;

  arma::field<arma::mat> F(3,1);
  F(0,0) = Uinit;
  F(1,0) = Vinit;
  F(2,0) = Laminit;

  unsigned int iter = 0;
  while( (err>tol) & (iter < maxiter)  ){
    iter++;

    for(unsigned int j = 0; j < m; j++) {
      if(j==0) {
        F(0,0).col(j) = X.col(j) + F(2,0).col(j);
      }
      else if(j < m-1) {
        F(0,0).col(j) = X.col(j) + F(2,0).col(j) - F(2,0).col(j-1);
      }
      else{
        F(0,0).col(j) = X.col(j) - F(2,0).col(j-1);
      }
    }
    ferrnew = 0;
    for(unsigned int k = 0; k < m-1; k++) {
      tmp = F(0,0).col(k) - F(0,0).col(k+1) - (1/nu)*F(2,0).col(k);
      tmpmax(0) = (1 - (sigma(k) / arma::norm(tmp,"fro")) );
      F(1,0).col(k) = arma::max(tmpmax)*tmp;
      ferrnew = ferrnew + w(k)*arma::norm(F(1,0).col(k),"fro");
      F(2,0).col(k) = F(2,0).col(k) + nu*(F(1,0).col(k) - F(0,0).col(k) + F(0,0).col(k+1));
    }
    ferrnew = arma::norm(X - F(0,0),"fro") + gamma*ferrnew;
    err = abs(ferrnew - ferrold) / abs(ferrold);
    ferrold = ferrnew;
  }
  return(F);

}
