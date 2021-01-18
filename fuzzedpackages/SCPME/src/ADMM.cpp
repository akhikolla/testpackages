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
//' \url{https://mgallow.github.io/SCPME/}.
//'
//' @param S pxp sample covariance matrix (denominator n).
//' @param A option to provide user-specified matrix for penalty term. This matrix must have p columns. Defaults to identity matrix.
//' @param B option to provide user-specified matrix for penalty term. This matrix must have p rows. Defaults to identity matrix.
//' @param C option to provide user-specified matrix for penalty term. This matrix must have nrow(A) rows and ncol(B) columns. Defaults to identity matrix.
//' @param initOmega initialization matrix for Omega
//' @param initZ initialization matrix for Z2
//' @param initY initialization matrix for Y
//' @param lam postive tuning parameter for elastic net penalty.
//' @param alpha elastic net mixing parameter contained in [0, 1]. \code{0 = ridge, 1 = lasso}. Alpha must be a single value (cross validation across alpha not supported).
//' @param rho initial step size for ADMM algorithm.
//' @param tau optional constant used to ensure positive definiteness in Q matrix in algorithm
//' @param mu factor for primal and residual norms in the ADMM algorithm. This will be used to adjust the step size \code{rho} after each iteration.
//' @param tau_rho factor in which to increase step size \code{rho}.
//' @param iter_rho step size \code{rho} will be updated every \code{iter.rho} steps
//' @param crit criterion for convergence (\code{ADMM} or \code{loglik}). If \code{crit = loglik} then iterations will stop when the relative change in log-likelihood is less than \code{tol.abs}. Default is \code{ADMM} and follows the procedure outlined in Boyd, et al.
//' @param tol_abs absolute convergence tolerance. Defaults to 1e-4.
//' @param tol_rel relative convergence tolerance. Defaults to 1e-4.
//' @param maxit maximum number of iterations. Defaults to 1e4.
//' 
//' @return returns list of returns which includes:
//' \item{Iterations}{number of iterations.}
//' \item{lam}{optimal tuning parameter.}
//' \item{Omega}{estimated penalized precision matrix.}
//' \item{Z2}{estimated Z matrix.}
//' \item{Y}{estimated Y matrix.}
//' \item{rho}{estimated rho.}
//' 
//' @references
//' \itemize{
//' \item Boyd, Stephen, Neal Parikh, Eric Chu, Borja Peleato, Jonathan Eckstein, and others. 2011. 'Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers.' \emph{Foundations and Trends in Machine Learning} 3 (1). Now Publishers, Inc.: 1-122. \url{https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf}
//' \item Hu, Yue, Chi, Eric C, amd Allen, Genevera I. 2016. 'ADMM Algorithmic Regularization Paths for Sparse Statistical Machine Learning.' \emph{Splitting Methods in Communication, Imaging, Science, and Engineering}. Springer: 433-459.
//' \item Molstad, Aaron J., and Adam J. Rothman. (2017). 'Shrinking Characteristics of Precision Matrix Estimators. \emph{Biometrika.}. \url{https://doi.org/10.1093/biomet/asy023}
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
List ADMMc(const arma::mat &S, const arma::mat &A, const arma::mat &B, const arma::mat &C, const arma::mat &initOmega, const arma::mat &initZ, const arma::mat &initY, const double lam, const double alpha = 1, const double tau = 10, double rho = 2, const double mu = 10, const double tau_rho = 2, const int iter_rho = 10, std::string crit = "ADMM", const double tol_abs = 1e-4, const double tol_rel = 1e-4, const int maxit = 1e4){
  
  // allocate memory
  bool criterion = true;
  int p = S.n_cols, iter = 0;
  double s, r, eps1, eps2, lik, lik2, sgn, logdet, sqrt;
  s = r = eps1 = eps2 = lik = lik2 = sgn = logdet = 0;
  arma::mat Z(initZ), Z2(initZ), Y(initY), Omega(initOmega), G, AOB, ACB, BZA, BYA;
  sqrt = std::sqrt(C.n_cols*C.n_rows);
  AOB = A*Omega*B;
  ACB = A.t()*C*B.t();

  
  // loop until convergence
  while (criterion && (iter < maxit)){
    
    // update values
    iter++;
    Z = Z2;
    
    // compute G (1)
    G = rho*A.t()*(AOB - Z + Y/rho)*B.t() - rho*ACB;
    
    // ridge equation (2)
    // gather eigen values (spectral decomposition)
    Omega = RIDGEc(S + (G + G.t())/2 - rho*tau*Omega, rho*tau);
    
    // penalty equation (3)
    // soft-thresholding
    AOB = A*Omega*B;
    Z2 = rho*(AOB - C) + Y;
    softmatrixc(Z2, lam*alpha);
    Z2 /= (rho + (1 - alpha)*lam);
    
    // update Y (4)
    Y += rho*(AOB - Z2 - C);
    
    // update dual and primal residuals
    BYA = B*Y.t()*A;
    BZA = B*(Z2 - Z).t()*A;
    s = arma::norm(rho*(BZA + BZA.t())/2, "fro");
    r = arma::norm(AOB - Z2 - C, "fro");
    
    // update dual and primal epsilons
    eps1 = sqrt*tol_abs + tol_rel*std::max(std::max(arma::norm(AOB, "fro"), arma::norm(Z2, "fro")), arma::norm(C, "fro"));
    eps2 = p*tol_abs + tol_rel*arma::norm((BYA + BYA.t())/2, "fro");
    
    // check rho update
    if (iter % iter_rho == 0){
      
      // increase step size
      if (r/eps1 > mu*s/eps2){
        rho *= tau_rho;
      }
      
      // decrease step size
      if (s/eps2 > mu*r/eps1){
        rho *= 1/tau_rho;
      }
    }
    
    // stopping criterion
    if (crit == "loglik"){
      
      // compute penalized loglik (close enough)
      arma::log_det(logdet, sgn, Omega);
      lik2 = (-p/2)*(arma::accu(Omega % S) - logdet + lam*((1 - alpha)/2*arma::accu(arma::square(AOB - C)) + alpha*arma::accu(arma::abs(AOB - C))));
      criterion = (std::abs((lik2 - lik)/lik) >= tol_abs);
      lik = lik2;
      
    } else {
      
      // ADMM criterion
      criterion = (r >= eps1 || s >= eps2);
      
    }
    
    // R_CheckUserInterrupt
    if (iter % 1000 == 0){
      R_CheckUserInterrupt();
    }
  }
  
  return List::create(Named("Iterations") = iter,
                      Named("lam") = lam,
                      Named("Omega") = Omega,
                      Named("Z") = Z2,
                      Named("Y") = Y,
                      Named("rho") = rho);
  
}

