// [[Rcpp::interfaces(r)]]
#include <sqp.h>   
 
//' @title
//' Auxiliary functions
//' 
//' @description
//' Auxiliary functions
//' 
//' @inherit qp_solver
//' 
//' @references \insertAllCited{}
//' 
//' @noRd
// [[Rcpp::export(.solvers_slacked_dense)]]      
Rcpp::List slacked_dense(          
    arma::vec                 x,        // Optimization variable (initial value)
    arma::mat                 Q,        // Quadratic Multiplier
    arma::mat                 C_eq,     // Equality  Constraint Multiplier
    arma::mat                 C_ineq,   // Inquality Constraint Multiplier
    arma::vec                 l,        // Linear Multiplier
    const arma::vec           &t_eq,    // Equality   Constraint RHS
    arma::vec                 t_ineq,   // Inequality Constraint RHS (upper bound)
    const double              &penalty, // Penalty parameter for slack variables
    const double              tol       = 1e-7,
    const unsigned            max_iter  = 500,
    int                       dim_eq    = -1,
    int                       dim_ineq  = -1,
    int                       dim_Q     = -1,
    const unsigned            solver    =  0,
    const bool                fast      = false,
    const bool                all_slack = false,
    const bool                debug     = false)
{
  return(sqp::solvers::slacked(
      x,
      Q,
      C_eq,
      C_ineq,
      l,
      t_eq,
      t_ineq,
      penalty,
      tol,
      max_iter,
      dim_eq,
      dim_ineq,
      dim_Q,
      solver,
      fast,
      all_slack,
      debug));
}
