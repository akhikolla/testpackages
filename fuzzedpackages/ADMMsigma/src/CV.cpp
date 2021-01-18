// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <progress.hpp>
#include "Sigma.h"
#include "soft.h"

using namespace Rcpp;




//' @title K fold (c++)
//' @description creates vector of shuffled indices.
//' @param n number of elements.
//' @param K number of folds.
//' @keywords internal
//'
arma::vec kfold(const int &n, const int &K){
  
  // create sequence 1:n
  arma::vec indices = arma::linspace<arma::vec>(1, n, n);
  
  // assign number fold
  for (int i = 0; i < n; i ++){
    indices[i] = i % K;
  }
  
  // shuffle indices
  indices = arma::shuffle(indices);
  
  return indices;
  
}



//--------------------------------------------------------------------------------------------




//' @title CV ADMM penalized precision matrix estimation (c++)
//' @description Cross validation function for ADMMsigma.
//'
//' @param X option to provide a nxp matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
//' @param S option to provide a pxp sample covariance matrix (denominator n). If argument is \code{NULL} and \code{X} is provided instead then \code{S} will be computed automatically.
//' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order.
//' @param alpha elastic net mixing parameter contained in [0, 1]. \code{0 = ridge, 1 = lasso}. If a vector of parameters is provided, they should be in increasing order.
//' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
//' @param path option to return the regularization path. This option should be used with extreme care if the dimension is large. If set to TRUE, cores will be set to 1 and errors and optimal tuning parameters will based on the full sample. Defaults to FALSE.
//' @param rho initial step size for ADMM algorithm.
//' @param mu factor for primal and residual norms in the ADMM algorithm. This will be used to adjust the step size \code{rho} after each iteration.
//' @param tau_inc factor in which to increase step size \code{rho}
//' @param tau_dec factor in which to decrease step size \code{rho}
//' @param crit criterion for convergence (\code{ADMM} or \code{loglik}). If \code{crit = loglik} then iterations will stop when the relative change in log-likelihood is less than \code{tol.abs}. Default is \code{ADMM} and follows the procedure outlined in Boyd, et al.
//' @param tol_rel relative convergence tolerance. Defaults to 1e-4.
//' @param maxit maximum number of iterations. Defaults to 1e4.
//' @param adjmaxit adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged (for each \code{alpha}). This option is intended to be paired with \code{warm} starts and allows for "one-step" estimators. Defaults to 1e4.
//' @param K specify the number of folds for cross validation.
//' @param crit_cv cross validation criterion (\code{loglik} \code{penloglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
//' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
//' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
//' 
//' @return list of returns includes:
//' \item{lam}{optimal tuning parameter.}
//' \item{alpha}{optimal tuning parameter.}
//' \item{path}{array containing the solution path. Solutions will be ordered in ascending alpha values for each lambda.}
//' \item{min.error}{minimum average cross validation error (cv_crit) for optimal parameters.}
//' \item{avg.error}{average cross validation error (cv_crit) across all folds.}
//' \item{cv.error}{cross validation errors (cv_crit).}
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
List CV_ADMMc(const arma::mat &X, const arma::mat &S, const arma::colvec &lam, const arma::colvec &alpha, bool diagonal = false, bool path = false, double rho = 2, const double mu = 10, const double tau_inc = 2, const double tau_dec = 2, std::string crit = "ADMM", const double tol_abs = 1e-4, const double tol_rel = 1e-4, int maxit = 1e4, int adjmaxit = 1e4, int K = 5, std::string crit_cv = "loglik", std::string start = "warm", std::string trace = "progress") {
  
  // initialization
  int n, p = S.n_cols, l = lam.n_rows, a = alpha.n_rows, initmaxit = maxit;
  double sgn = 0, logdet = 0, initrho = rho, alpha_, lam_;
  arma::mat X_train, X_test, S_train(S), S_test(S), Omega, initOmega, initZ, initY;
  arma::mat CV_error, zerosla(l, a, arma::fill::zeros);
  arma::uvec index, index_; arma::vec folds; arma::rowvec X_bar;
  arma::cube CV_errors(l, a, K, arma::fill::zeros), Path;
  Progress progress(l*a*K, trace == "progress");
  
  // no need to create folds if K = 1
  if (K == 1){
    
    // set sample size
    n = S.n_rows;
    
    // initialize Path, if necessary
    if (path){
      Path = arma::zeros<arma::cube>(p, p, l*a);
    }
    
  } else {
    
    // designate folds and shuffle -- ensures randomized folds
    n = X.n_rows;
    folds = kfold(n, K);
    
  }
  
  // parse data into folds and perform CV
  for (int k = 0; k < K; k++){
    
    // re-initialize values for each fold
    CV_error = zerosla; maxit = initmaxit; rho = initrho;
    initOmega = initZ = arma::diagmat(1/arma::diagvec(S));
    initY = arma::zeros<arma::mat>(p, p);
      
    if (K > 1) {
      
      // separate into training and testing data
      index = arma::find(folds != k);
      index_ = arma::find(folds == k);
      
      // training set
      X_train = X.rows(index);
      X_bar = arma::mean(X_train, 0);
      X_train -= arma::ones<arma::colvec>(X_train.n_rows)*X_bar;
      
      // validation set
      X_test = X.rows(index_);
      X_test -= arma::ones<arma::colvec>(X_test.n_rows)*X_bar;
      n = X_test.n_rows;
      
      // sample covariances
      S_train = arma::cov(X_train, 1);
      S_test = arma::cov(X_test, 1);
      
    }
    
    // loop over all tuning parameters
    for (int i = 0; i < l; i++){
      for (int j = 0; j < a; j++){
        
        // set temporary tuning parameters
        lam_ = lam[i];
        alpha_ = alpha[j];
        
        // compute the ridge-penalized likelihood precision matrix estimator at the ith value in lam:
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
        CV_error(i, j) = (n/2)*(arma::accu(Omega % S_test) - logdet);
        
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
        
        // save estimate if path = TRUE
        if (path){
          Path.slice(j + i*a) = as<arma::mat>(ADMM["Z"]);
        }
        
        // update progress bar
        if (trace == "progress"){
          progress.increment();
          
          // if not quiet, then print progress lambda
        } else if (trace == "print"){
          Rcout << "Finished lam = " << lam[i] << " in fold " << k << "\n";
        }
      }
    }
    
    // if not quiet, then print progress fold
    if (trace == "print"){
      Rcout << "Finished fold" << k << "\n";
    }
    
    // append CV errors
    CV_errors.slice(k) = CV_error;
    
  }

  
  // determine optimal tuning parameters
  arma::mat AVG_error = arma::mean(CV_errors, 2);
  double error = AVG_error.min();
  arma::uword ind = AVG_error.index_min();
  int lam_ind = ind % AVG_error.n_rows;
  int alpha_ind = std::floor(ind/AVG_error.n_rows);
  double best_lam = lam[lam_ind];
  double best_alpha = alpha[alpha_ind];
  
  
  // return list of coefficients
  return List::create(Named("lam") = best_lam,
                      Named("alpha") = best_alpha,
                      Named("path") = Path,
                      Named("min.error") = error,
                      Named("avg.error") = AVG_error,
                      Named("cv.error") = CV_errors);
}



////-----------------------------------------------------




//' @title CV ridge penalized precision matrix estimation (c++)
//' @description Cross validation function for RIDGEsigma.
//' 
//' @param X option to provide a nxp matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
//' @param S option to provide a pxp sample covariance matrix (denominator n). If argument is \code{NULL} and \code{X} is provided instead then \code{S} will be computed automatically.
//' @param lam positive tuning parameters for ridge penalty. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{10^seq(-5, 5, 0.5)}.
//' @param path option to return the regularization path. This option should be used with extreme care if the dimension is large. If set to TRUE, cores will be set to 1 and errors and optimal tuning parameters will based on the full sample. Defaults to FALSE.
//' @param K specify the number of folds for cross validation.
//' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
//' 
//' @return list of returns includes:
//' \item{lam}{optimal tuning parameter.}
//' \item{path}{array containing the solution path. Solutions are ordered dense to sparse.}
//' \item{min.error}{minimum average cross validation error for optimal parameters.}
//' \item{avg.error}{average cross validation error across all folds.}
//' \item{cv.error}{cross validation errors (negative validation likelihood).}
//'
//' @keywords internal
//'
// [[Rcpp::export]]
List CV_RIDGEc(const arma::mat &X, const arma::mat &S, const arma::colvec &lam, bool path = false, int K = 3, std::string trace = "none") {
  
  // initialization
  int n, p = S.n_cols, l = lam.n_rows;
  double sgn = 0, logdet = 0, lam_;
  arma::mat X_train, X_test, S_train(S), S_test(S);
  arma::mat Omega, CV_errors(l, K, arma::fill::zeros), zeros(p, p, arma::fill::zeros);
  arma::uvec index, index_; arma::vec folds; arma::cube Path;
  arma::colvec CV_error, zerosl(l, arma::fill::zeros); arma::rowvec X_bar;
  Progress progress(l*K, trace == "progress");
  
  // no need to create folds if K = 1
  if (K == 1){
    
    // set sample size
    n = S.n_rows;
    
    // initialize Path, if necessary
    if (path){
      Path = arma::zeros<arma::cube>(p, p, l);
    }
    
  } else {
    
    // designate folds and shuffle -- ensures randomized folds
    n = X.n_rows;
    folds = kfold(n, K);
    
  }
  
  // parse data into folds and perform CV
  for (int k = 0; k < K; k++){
    
    // re-initialize for each fold
    CV_error = zerosl;
      
    if (K > 1) {
    
    // separate into training and testing data
    index = arma::find(folds != k);
    index_ = arma::find(folds == k);
    
    // training set
    X_train = X.rows(index);
    X_bar = arma::mean(X_train, 0);
    X_train -= arma::ones<arma::colvec>(X_train.n_rows)*X_bar;
    
    // validation set
    X_test = X.rows(index_);
    X_test -= arma::ones<arma::colvec>(X_test.n_rows)*X_bar;
    n = X_test.n_rows;
    
    // sample covariances
    S_train = arma::cov(X_train, 1);
    S_test = arma::cov(X_test, 1);
    
    }
    
    // loop over all tuning parameters
    for (int i = 0; i < l; i++){
      
      // set temporary tuning parameters
      lam_ = lam[i];
      
      // compute the ridge-penalized likelihood precision matrix estimator at the ith value in lam:
      Omega = RIDGEc(S_train, lam_);
      
      // compute the observed negative validation loglikelihood (close enough)
      arma::log_det(logdet, sgn, Omega);
      CV_error[i] = (n/2)*(arma::accu(Omega % S_test) - logdet);
      
      // save estimate if path = TRUE
      if (path){
        Path.slice(i) = Omega;
      }
      
      // update progress bar
      if (trace == "progress"){
        progress.increment();
      
      // if not quiet, then print progress lambda
      } else if (trace == "print"){
        Rcout << "Finished lam = " << lam[i] << " in fold " << k << "\n";
      }
    }
    
    if (trace == "print"){
      Rcout << "Finished fold" << k << "\n";
    }
    
    // append CV errors
    CV_errors.col(k) = CV_error;
    
  }
  
  // determine optimal tuning parameters
  arma::colvec AVG_error = arma::mean(CV_errors, 1);
  double error = AVG_error.min();
  arma::uword ind = AVG_error.index_min();
  int lam_ind = ind % AVG_error.n_rows;
  double best_lam = lam[lam_ind];
  
  // return list of coefficients
  return List::create(Named("lam") = best_lam,
                      Named("path") = Path,
                      Named("min.error") = error,
                      Named("avg.error") = AVG_error,
                      Named("cv.error") = CV_errors);
}
