// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <progress.hpp>
#include "ADMM.h"
#include "soft.h"

using namespace Rcpp;




//' @title CV (no folds) ADMM penalized precision matrix estimation (c++)
//' @description Cross validation (no folds) function for shrink. This function is to be used with CVP_ADMM.
//'
//' @param X_train nxp training data matrix.
//' @param X_valid (n - q)xp validation data matrix matrix.
//' @param Y_train nxr training response matrix.
//' @param Y_valid (n - q)xr validation response matrix.
//' @param A option to provide user-specified matrix for penalty term. This matrix must have p columns. Defaults to identity matrix.
//' @param B option to provide user-specified matrix for penalty term. This matrix must have p rows. Defaults to identity matrix.
//' @param C option to provide user-specified matrix for penalty term. This matrix must have nrow(A) rows and ncol(B) columns. Defaults to identity matrix.
//' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order.
//' @param alpha elastic net mixing parameter contained in [0, 1]. \code{0 = ridge, 1 = lasso}. Alpha must be a single value (cross validation across alpha not supported).
//' @param tau optional constant used to ensure positive definiteness in Q matrix in algorithm
//' @param rho initial step size for ADMM algorithm.
//' @param mu factor for primal and residual norms in the ADMM algorithm. This will be used to adjust the step size \code{rho} after each iteration.
//' @param tau_rho factor in which to increase/decrease step size \code{rho}
//' @param iter_rho step size \code{rho} will be updated every \code{iter.rho} steps
//' @param crit criterion for convergence (\code{ADMM} or \code{loglik}). If \code{crit = loglik} then iterations will stop when the relative change in log-likelihood is less than \code{tol.abs}. Default is \code{ADMM} and follows the procedure outlined in Boyd, et al.
//' @param tol_abs absolute convergence tolerance. Defaults to 1e-4.
//' @param tol_rel relative convergence tolerance. Defaults to 1e-4.
//' @param maxit maximum number of iterations. Defaults to 1e4.
//' @param adjmaxit adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged. This option is intended to be paired with \code{warm} starts and allows for "one-step" estimators. Defaults to 1e4.
//' @param crit_cv cross validation criterion (\code{MSE}, \code{loglik}, \code{penloglik} \code{AIC}, or \code{BIC}). Defaults to \code{MSE}.
//' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
//' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
//' 
//' @return cross validation errors (cv_crit)
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat CVP_ADMMc(const arma::mat &X_train, const arma::mat &X_valid, const arma::mat &Y_train, const arma::mat &Y_valid, const arma::mat &A, const arma::mat &B, const arma::mat &C, const arma::colvec &lam, const double alpha = 1, const double tau = 10, double rho = 2, const double mu = 10, const double tau_rho = 2, const int iter_rho = 10, std::string crit = "ADMM", const double tol_abs = 1e-4, const double tol_rel = 1e-4, int maxit = 1e4, int adjmaxit = 1e4, std::string crit_cv = "MSE", std::string start = "warm", std::string trace = "progress") {
  
  // initialization
  int n = X_valid.n_rows, r = Y_valid.n_cols, l = lam.n_rows;
  double sgn = 0, logdet = 0, lam_;
  arma::mat S_train, S_valid, Omega, initOmega, initZ, initY; arma::colvec nzeros;
  arma::mat CV_error(l, 1, arma::fill::zeros);
  Progress progress(l, trace == "progress");
  
  // calculate sample covariances
  S_train = arma::cov(X_train, 1);
  S_valid = arma::cov(X_valid, 1);
  
  // initial estimates
  initOmega = arma::diagmat(1/arma::diagvec(S_train));
  initZ = A*initOmega*B - C;
  initY = arma::zeros<arma::mat>(C.n_rows, C.n_cols);
  
  // loop over all tuning parameters
  for (int i = 0; i < l; i++){
    
    // set temporary tuning parameters
    lam_ = lam[i];
    
    // compute the penalized likelihood precision matrix estimator at the ith value in lam:
    List ADMM = ADMMc(S_train, A, B, C, initOmega, initZ, initY, lam_, alpha, tau, rho, mu, tau_rho, iter_rho, crit, tol_abs, tol_rel, maxit);
    Omega = as<arma::mat>(ADMM["Omega"]);
    
    if (start == "warm"){
      
      // option to save initial values for warm starts
      initOmega = as<arma::mat>(ADMM["Omega"]);
      initZ = as<arma::mat>(ADMM["Z"]);
      initY = as<arma::mat>(ADMM["Y"]);
      rho = as<double>(ADMM["rho"]);
      maxit = adjmaxit;
      
    }
    
    // criterion MSE
    if (crit_cv == "MSE"){
      CV_error[i] = std::pow(arma::norm(Y_valid - X_valid*Omega*arma::cov(X_valid, Y_valid, 1), "fro"), 2)/(n*r);
      
    // criterion loglik
    } else if (crit_cv == "loglik"){
      arma::log_det(logdet, sgn, Omega);
      CV_error[i] = (n/2)*(arma::accu(Omega % S_valid) - logdet);
      
    // critertion penloglik
    } else if (crit_cv == "penloglik"){
      arma::log_det(logdet, sgn, Omega);
      CV_error[i] = (n/2)*(arma::accu(Omega % S_valid) - logdet) + lam_*((1 - alpha)/2*arma::accu(arma::square(A*Omega*B - C)) + alpha*arma::accu(arma::abs(A*Omega*B - C)));
      
      
    // criterion AIC
    } else if (crit_cv == "AIC"){
      arma::log_det(logdet, sgn, Omega);
      CV_error[i] = (n/2)*(arma::accu(Omega % S_valid) - logdet) + numzeros(Omega);
      
    // criterion BIC
    } else {
      arma::log_det(logdet, sgn, Omega);
      CV_error[i] = (n/2)*(arma::accu(Omega % S_valid) - logdet) + numzeros(Omega)*std::log(n)/2;
    }
    
    // update progress bar
    if (trace == "progress"){
      progress.increment();
      
    // if not quiet, then print progress lambda
    } else if (trace == "print"){
      Rcout << "Finished lam = " << lam[i] << "\n";
    }
  }
  
  // return CV errors
  return(CV_error);
}


