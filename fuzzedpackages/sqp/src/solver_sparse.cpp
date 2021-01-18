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
// [[Rcpp::export(.solvers_slacked_sparse)]] 
Rcpp::List slacked_sparse(
    arma::vec                 x,        // Optimization variable (initial value)
    arma::sp_mat              Q,        // Quadratic Multiplier
    arma::sp_mat              C_eq,     // Equality  Constraint Multiplier
    arma::sp_mat              C_ineq,   // Inquality Constraint Multiplier
    arma::vec                 l,        // Linear Multiplier
    const arma::vec           &t_eq,    // Equality   Constraint RHS
    arma::vec                 t_ineq,   // Inequality Constraint RHS (upper bound)
    const double              &penalty, // Penalty parameter for slack variables
    const double              tol       = 1e-7,
    const unsigned            max_iter  = 500,
    int                       dim_eq    = -1,
    int                       dim_ineq  = -1,
    int                       dim_Q     = -1,
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
      all_slack,
      debug));
}
