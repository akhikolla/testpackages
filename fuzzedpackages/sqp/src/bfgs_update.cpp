// [[Rcpp::interfaces(r)]] 
#include <sqp.h>   

//' @title
//' (Damped) BFGS Hessian approximation
//'  
//' @description
//' \code{BFGS} update for appromation of the Hessian matrix
//' \insertCite{@cf. @broyden70; @fletcher70; @goldfarb70; @shanno70}{sqp}
//' in its damped version proposed by \insertCite{powell78;textual}{sqp}.
//' The approximation is based on first-order information  
//' (parameter values & gradients) only.
//' 
//' @param hessian
//' \strong{Dense matrix}
//' of size \eqn{N \times N}{N x N}:\cr
//' Current approximation of the Hessian matrix, which is
//' updated by reference.\cr
//' Needs to be symmetric positive definite.\cr 
//' A common starting point for the \code{BFGS} algorithm is the identity matrix.
//'                                   
//' @param  old_y,new_y,old_gradient,new_gradient
//' \strong{Numeric vectors} 
//' of size \code{N}:\cr
//' parameters \bold{old_y},\bold{new_y} and\cr
//' corresponding gradients \bold{old_gradient},\bold{new_gradient}
//' from previous and current iteration.
//' 
//' @param constraint_adjustment
//' \strong{Boolean}:\cr
//' Whether to enforce positive definiteness\cr
//' (mainly for constrained optimization).
//' 
//' @return 
//' Nothing. Argument 'hessian' is updated by reference.
//' 
//' @references \insertAllCited{}
//' 
//' @export
// [[Rcpp::export]] 
void bfgs_update(arma::mat &hessian,
                 arma::vec &old_y,
                 arma::vec &new_y,
                 arma::vec &old_gradient,
                 arma::vec &new_gradient,
                 const bool constraint_adjustment = true)
{
  sqp::bfgs::bfgs_update(hessian,
                        old_y,
                        new_y,
                        old_gradient,
                        new_gradient,
                        constraint_adjustment);
}
