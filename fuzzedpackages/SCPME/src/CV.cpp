// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <progress.hpp>
#include "ADMM.h"
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
//' @param Y option to provide nxr response matrix. Each row corresponds to a single response and each column contains n response of a single feature/response.
//' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order.
//' @param alpha elastic net mixing parameter contained in [0, 1]. \code{0 = ridge, 1 = lasso}. Alpha must be a single value (cross validation across alpha not supported).
//' @param A option to provide user-specified matrix for penalty term. This matrix must have p columns. Defaults to identity matrix.
//' @param B option to provide user-specified matrix for penalty term. This matrix must have p rows. Defaults to identity matrix.
//' @param C option to provide user-specified matrix for penalty term. This matrix must have nrow(A) rows and ncol(B) columns. Defaults to identity matrix.
//' @param path option to return the regularization path. This option should be used with extreme care if the dimension is large. If set to TRUE, cores will be set to 1 and errors and optimal tuning parameters will based on the full sample. Defaults to FALSE.
//' @param tau optional constant used to ensure positive definiteness in Q matrix in algorithm
//' @param rho initial step size for ADMM algorithm.
//' @param mu factor for primal and residual norms in the ADMM algorithm. This will be used to adjust the step size \code{rho} after each iteration.
//' @param tau_rho factor in which to increase/decrease step size \code{rho}
//' @param iter_rho step size \code{rho} will be updated every \code{iter.rho} steps
//' @param crit criterion for convergence (\code{ADMM} or \code{loglik}). If \code{crit = loglik} then iterations will stop when the relative change in log-likelihood is less than \code{tol.abs}. Default is \code{ADMM} and follows the procedure outlined in Boyd, et al.
//' @param tol_rel relative convergence tolerance. Defaults to 1e-4.
//' @param maxit maximum number of iterations. Defaults to 1e4.
//' @param adjmaxit adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged. This option is intended to be paired with \code{warm} starts and allows for "one-step" estimators. Defaults to 1e4.
//' @param K specify the number of folds for cross validation.
//' @param crit_cv cross validation criterion (\code{loglik}, \code{penloglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
//' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
//' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
//' 
//' @return list of returns includes:
//' \item{lam}{optimal tuning parameter.}
//' \item{path}{array containing the solution path. Solutions will be ordered in ascending lambda values.}
//' \item{min.error}{minimum average cross validation error (cv_crit) for optimal parameters.}
//' \item{avg.error}{average cross validation error (cv_crit) across all folds.}
//' \item{cv.error}{cross validation errors (cv_crit).}
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
List CV_ADMMc(const arma::mat &X, const arma::mat &S, const arma::mat &Y, const arma::mat &A, const arma::mat &B, const arma::mat &C, const arma::colvec &lam, const double alpha = 1, bool path = false, double tau = 10, double rho = 2, const double mu = 10, const double tau_rho = 2, const int iter_rho = 10, std::string crit = "ADMM", const double tol_abs = 1e-4, const double tol_rel = 1e-4, int maxit = 1e4, int adjmaxit = 1e4, int K = 5, std::string crit_cv = "MSE", std::string start = "warm", std::string trace = "progress") {
  
  // initialization
  int n, r = Y.n_cols, p = S.n_cols, l = lam.n_rows, initmaxit = maxit;
  double sgn = 0, logdet = 0, initrho = rho, lam_;
  arma::mat X_train, X_valid, Y_train, Y_valid, S_train(S), S_valid(S), Omega, initOmega, initZ, initY;
  arma::mat CV_errors(l, K, arma::fill::zeros);
  arma::colvec CV_error, zerosl(l, arma::fill::zeros); arma::rowvec X_bar;
  arma::uvec index, index_; arma::vec folds; arma::cube Path;
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
    
    // re-initialize values for each fold
    CV_error = zerosl; maxit = initmaxit; rho = initrho;
    initOmega = arma::diagmat(1/arma::diagvec(S));
    initZ = A*initOmega*B - C;
    initY = arma::zeros<arma::mat>(C.n_rows, C.n_cols);
      
    if (K > 1) {
      
      // separate into training and testing data
      index = arma::find(folds != k);
      index_ = arma::find(folds == k);
      
      // training set
      X_train = X.rows(index);
      X_bar = arma::mean(X_train, 0);
      X_train -= arma::ones<arma::colvec>(X_train.n_rows)*X_bar;
      
      // validation set
      X_valid = X.rows(index_);
      X_valid -= arma::ones<arma::colvec>(X_valid.n_rows)*X_bar;
      n = X_valid.n_rows;
      
      // sample covariances
      S_train = arma::cov(X_train, 1);
      S_valid = arma::cov(X_valid, 1);
      
      // training/validation for Y, if necessary
      if (crit_cv == "MSE"){
        
        Y_train = Y.rows(index);
        Y_valid = Y.rows(index_);
        
      }
      
    }
    
    // loop over all tuning parameters
    for (int i = 0; i < l; i++){
        
        // set temporary tuning parameters
        lam_ = lam[i];
        
        // compute the ridge-penalized likelihood precision matrix estimator at the ith value in lam:
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
          CV_error[i] = arma::accu(arma::square(Y_valid - X_valid*Omega*arma::cov(X_valid, Y_valid, 1)))/(n*r);
        
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
        
        // save estimate if path = TRUE
        if (path){
          Path.slice(i) = as<arma::mat>(ADMM["Omega"]);
        }
        
        // update progress bar
        if (trace == "progress"){
          progress.increment();
          
          // if not quiet, then print progress lambda
        } else if (trace == "print"){
          Rcout << "Finished lam = " << lam[i] << " in fold " << k << "\n";
        }
      }
    
    // if not quiet, then print progress fold
    if (trace == "print"){
      Rcout << "Finished fold" << k << "\n";
    }
    
    // append CV errors
    CV_errors.col(k) = CV_error;
    
  }

  
  // determine optimal tuning parameters
  arma::mat AVG_error = arma::mean(CV_errors, 1);
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


