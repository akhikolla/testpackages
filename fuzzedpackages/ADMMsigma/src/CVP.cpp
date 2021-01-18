// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <progress.hpp>
#include "Sigma.h"
#include "soft.h"

using namespace Rcpp;




//' @title CV (no folds) ADMM penalized precision matrix estimation (c++)
//' @description Cross validation (no folds) function for ADMMsigma. This function is to be used with CVP_ADMM.
//'
//' @param n sample size for X_valid (used to calculate crit_cv)
//' @param S_train pxp sample covariance matrix for training data (denominator n).
//' @param S_valid pxp sample covariance matrix for validation data (denominator n).
//' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order.
//' @param alpha elastic net mixing parameter contained in [0, 1]. \code{0 = ridge, 1 = lasso}. If a vector of parameters is provided, they should be in increasing order.
//' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
//' @param rho initial step size for ADMM algorithm.
//' @param mu factor for primal and residual norms in the ADMM algorithm. This will be used to adjust the step size \code{rho} after each iteration.
//' @param tau_inc factor in which to increase step size \code{rho}
//' @param tau_dec factor in which to decrease step size \code{rho}
//' @param crit criterion for convergence (\code{ADMM} or \code{loglik}). If \code{crit = loglik} then iterations will stop when the relative change in log-likelihood is less than \code{tol.abs}. Default is \code{ADMM} and follows the procedure outlined in Boyd, et al.
//' @param tol_abs absolute convergence tolerance. Defaults to 1e-4.
//' @param tol_rel relative convergence tolerance. Defaults to 1e-4.
//' @param maxit maximum number of iterations. Defaults to 1e4.
//' @param adjmaxit adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged (for each \code{alpha}). This option is intended to be paired with \code{warm} starts and allows for "one-step" estimators. Defaults to 1e4.
//' @param crit_cv cross validation criterion (\code{loglik}, \code{penloglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
//' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
//' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
//' 
//' @return cross validation errors (cv_crit)
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat CVP_ADMMc(const int n, const arma::mat &S_train, const arma::mat &S_valid, const arma::colvec &lam, const arma::colvec &alpha, bool diagonal = false, double rho = 2, const double mu = 10, const double tau_inc = 2, const double tau_dec = 2, std::string crit = "ADMM", const double tol_abs = 1e-4, const double tol_rel = 1e-4, int maxit = 1e4, int adjmaxit = 1e4, std::string crit_cv = "loglik", std::string start = "warm", std::string trace = "progress") {
  
  // initialization
  int p = S_train.n_rows, l = lam.n_rows, a = alpha.n_rows;
  double sgn = 0, logdet = 0, alpha_, lam_;
  arma::mat Omega, initOmega, initZ, initY; arma::colvec nzeros;
  initOmega = initZ = arma::diagmat(1/arma::diagvec(S_train));
  initY = arma::zeros<arma::mat>(p, p);
  arma::mat CV_error(l, a, arma::fill::zeros);
  Progress progress(l*a, trace == "progress");
  
  
  // loop over all tuning parameters
  for (int i = 0; i < l; i++){
    for (int j = 0; j < a; j++){
      
      // set temporary tuning parameters
      lam_ = lam[i];
      alpha_ = alpha[j];
      
      // compute the penalized likelihood precision matrix estimator at the ith value in lam:
      List ADMM = ADMMc(S_train, initOmega, initZ, initY, lam_, alpha_, diagonal, rho, mu, tau_inc, tau_dec, crit, tol_abs, tol_rel, maxit);
      Omega = as<arma::mat>(ADMM["Omega"]);

      if (start == "warm"){
        
        // option to save initial values for warm starts
        initOmega = as<arma::mat>(ADMM["Omega"]);
        initZ = as<arma::mat>(ADMM["Z"]);
        initY = as<arma::mat>(ADMM["Y"]);
        rho = as<double>(ADMM["rho"]);
        maxit = adjmaxit;
        
      }
      
      // compute the observed negative validation loglikelihood (close enough)
      arma::log_det(logdet, sgn, Omega);
      CV_error(i, j) = (n/2)*(arma::accu(Omega % S_valid) - logdet);
      
      // update for crit_cv, if necessary
      if (crit_cv == "penloglik"){
        CV_error(i, j) += lam_*((1 - alpha_)/2*arma::accu(arma::square(Omega)) + alpha_*arma::accu(arma::abs(Omega)));
      }
      if (crit_cv == "AIC"){
        CV_error(i, j) += numzeros(Omega);
      }
      if (crit_cv == "BIC"){
        CV_error(i, j) += numzeros(Omega)*std::log(n)/2;
      }
      
      // update progress bar
      if (trace == "progress"){
        progress.increment();
      
      // if not quiet, then print progress lambda
      } else if (trace == "print"){
        Rcout << "Finished lam = " << lam[i] << "\n";
      }
    }
  }
  
  // return CV errors
  return(CV_error);
}





//-------------------------------------------------------------------------------------





//' @title CV (no folds) RIDGE penalized precision matrix estimation (c++)
//' @description Cross validation (no folds) function for RIDGEsigma. This function is to be used with CVP_RIDGE.
//'
//' @param n sample size for X_valid (used to calculate CV_error)
//' @param S_train pxp sample covariance matrix for training data (denominator n).
//' @param S_valid pxp sample covariance matrix for validation data (denominator n).
//' @param lam positive tuning parameters for ridge penalty. If a vector of parameters is provided, they should be in increasing order.
//' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
//' 
//' @return cross validation errors (negative validation likelihood)
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat CVP_RIDGEc(const int n, const arma::mat &S_train, const arma::mat &S_valid, const arma::colvec &lam, std::string trace = "none") {
  
  // initialization
  int p = S_train.n_rows, l = lam.n_rows;
  double sgn = 0, logdet = 0, lam_;
  arma::mat Omega(p, p, arma::fill::zeros), CV_error(l, 1, arma::fill::zeros);
  Progress progress(l, trace == "progress");
  
  
  // loop over all tuning parameters
  for (int i = 0; i < l; i++){
    
    // set temporary tuning parameters
    lam_ = lam[i];
    
    // compute the ridge-penalized likelihood precision matrix estimator at the ith value in lam:
    Omega = RIDGEc(S_train, lam_);
    
    // compute the observed negative validation loglikelihood (close enough)
    arma::log_det(logdet, sgn, Omega);
    CV_error[i] = (n/2)*(arma::accu(Omega % S_valid) - logdet);
    
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
