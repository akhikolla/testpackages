// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "soft.h"

using namespace Rcpp;



//' @title Ridge-penalized precision matrix estimation (c++)
//' @description Ridge penalized matrix estimation via closed-form solution. Augmented from Adam Rothman's STAT 8931 code.
//'
//' @param S sample covariance matrix (denominator n).
//' @param lam tuning parameter for ridge penalty.
//' 
//' @return estimated Omega
//' 
//' @export
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat RIDGEc(const arma::mat &S, double lam){

  // gather eigen values of S (spectral decomposition)
  arma::mat V;
  arma::colvec Q;
  eig_sym(Q, V, S);

  // augment eigen values for omega hat
  arma::mat Q2 = (-Q + arma::sqrt(arma::square(Q) + 4*lam))/(2*lam);

  // compute omega hat for lambda (zero gradient equation)
  arma::mat omega = V*arma::diagmat(Q2)*V.t();

  return(omega);

}



////-----------------------------------------------------



//' @title Penalized precision matrix estimation via ADMM (c++)
//' 
//' @description Penalized precision matrix estimation using the ADMM algorithm
//' 
//' @details For details on the implementation of 'ADMMsigma', see the vignette
//' \url{https://mgallow.github.io/ADMMsigma/}.
//'
//' @param S pxp sample covariance matrix (denominator n).
//' @param initOmega initialization matrix for Omega
//' @param initZ initialization matrix for Z
//' @param initY initialization matrix for Y
//' @param lam postive tuning parameter for elastic net penalty.
//' @param alpha elastic net mixing parameter contained in [0, 1]. \code{0 = ridge, 1 = lasso}. Defaults to alpha = 1.
//' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
//' @param rho initial step size for ADMM algorithm.
//' @param mu factor for primal and residual norms in the ADMM algorithm. This will be used to adjust the step size \code{rho} after each iteration.
//' @param tau_inc factor in which to increase step size \code{rho}.
//' @param tau_dec factor in which to decrease step size \code{rho}.
//' @param crit criterion for convergence (\code{ADMM} or \code{loglik}). If \code{crit = loglik} then iterations will stop when the relative change in log-likelihood is less than \code{tol.abs}. Default is \code{ADMM} and follows the procedure outlined in Boyd, et al.
//' @param tol_abs absolute convergence tolerance. Defaults to 1e-4.
//' @param tol_rel relative convergence tolerance. Defaults to 1e-4.
//' @param maxit maximum number of iterations. Defaults to 1e4.
//' 
//' @return returns list of returns which includes:
//' \item{Iterations}{number of iterations.}
//' \item{lam}{optimal tuning parameters.}
//' \item{alpha}{optimal tuning parameter.}
//' \item{Omega}{estimated penalized precision matrix.}
//' \item{Z2}{estimated Z matrix.}
//' \item{Y}{estimated Y matrix.}
//' \item{rho}{estimated rho.}
//' 
//' @references
//' \itemize{
//' \item Boyd, Stephen, Neal Parikh, Eric Chu, Borja Peleato, Jonathan Eckstein, and others. 2011. 'Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers.' \emph{Foundations and Trends in Machine Learning} 3 (1). Now Publishers, Inc.: 1-122. \url{https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf}
//' \item Hu, Yue, Chi, Eric C, amd Allen, Genevera I. 2016. 'ADMM Algorithmic Regularization Paths for Sparse Statistical Machine Learning.' \emph{Splitting Methods in Communication, Imaging, Science, and Engineering}. Springer: 433-459.
//' \item Zou, Hui and Hastie, Trevor. 2005. "Regularization and Variable Selection via the Elastic Net." \emph{Journal of the Royal Statistial Society: Series B (Statistical Methodology)} 67 (2). Wiley Online Library: 301-320.
//' \item Rothman, Adam. 2017. 'STAT 8931 notes on an algorithm to compute the Lasso-penalized Gaussian likelihood precision matrix estimator.'
//' }
//' 
//' @author Matt Galloway \email{gall0441@@umn.edu}
//' 
//' @export
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
List ADMMc(const arma::mat &S, const arma::mat &initOmega, const arma::mat &initZ, const arma::mat &initY, const double lam, const double alpha = 1, bool diagonal = false, double rho = 2, const double mu = 10, const double tau_inc = 2, const double tau_dec = 2, std::string crit = "ADMM", const double tol_abs = 1e-4, const double tol_rel = 1e-4, const int maxit = 1e4){
  
  // allocate memory
  bool criterion = true;
  int p = S.n_cols, iter = 0;
  double s, r, eps1, eps2, lik, lik2, sgn, logdet;
  s = r = eps1 = eps2 = lik = lik2 = sgn = logdet = 0;
  arma::mat Z2(initZ), Z(initZ), Y(initY), Omega(initOmega), C(p, p, arma::fill::ones), Tau, Taum;

  
  // option to penalize diagonal elements
  if (!diagonal){
    C -= arma::diagmat(C);
  }
  
  // save values
  Tau = lam*alpha*C;
  Taum = lam*C - Tau;
  
  // loop until convergence
  while (criterion && (iter < maxit)){
    
    // update values
    iter++;
    Z = Z2;
    
    // ridge equation (1)
    // gather eigen values (spectral decomposition)
    Omega = RIDGEc(S + Y - rho*Z, rho);
    
    // penalty equation (2)
    // soft-thresholding
    Z2 = Y + rho*Omega;
    softmatrixc(Z2, Tau);
    Z2 /= (Taum + rho);
    
    // update Y (3)
    Y += rho*(Omega - Z2);
    
    // calculate new rho
    s = arma::norm(rho*(Z2 - Z), "fro");
    r = arma::norm(Omega - Z2, "fro");
    if (r > mu*s){
      rho *= tau_inc;
    }
    if (s > mu*r){
      rho *= 1/tau_dec;
    }
    
    // stopping criterion
    if (crit == "loglik"){
      
      // compute penalized loglik (close enough)
      arma::log_det(logdet, sgn, Omega);
      lik2 = (-p/2)*(arma::accu(Omega % S) - logdet + lam*((1 - alpha)/2*arma::accu(arma::square(C % Omega)) + alpha*arma::accu(C % arma::abs(Omega))));
      criterion = (std::abs((lik2 - lik)/lik) >= tol_abs);
      lik = lik2;
      
    } else {
      
      // ADMM criterion
      eps1 = p*tol_abs + tol_rel*std::max(arma::norm(Omega, "fro"), arma::norm(Z2, "fro"));
      eps2 = p*tol_abs + tol_rel*arma::norm(Y, "fro");
      criterion = (r >= eps1 || s >= eps2);
      
    }
    
    // R_CheckUserInterrupt
    if (iter % 1000 == 0){
      R_CheckUserInterrupt();
    }
  }
  
  return List::create(Named("Iterations") = iter,
                      Named("lam") = lam,
                      Named("alpha") = alpha,
                      Named("Omega") = Omega,
                      Named("Z") = Z2,
                      Named("Y") = Y,
                      Named("rho") = rho);
  
}

