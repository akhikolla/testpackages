#ifndef __sqp_solvers_quadratic_slack_included__        // if include guard for 'solvers/quadratic_slack.h' is undefined
#define __sqp_solvers_quadratic_slack_included__        // define include guard for 'solvers/quadratic_slack.h'  

#include "sqp.h"

namespace sqp {
namespace solvers{



inline sqp::solvers::result slacked(
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
  // Rcpp::Rcout << "dim_eq: " << dim_eq << "\n";
  // Rcpp::Rcout << "dim_ineq: " << dim_ineq << "\n";
  if(dim_eq < 0)
    dim_eq = C_eq.n_rows;
  
  if(dim_ineq < 0)
    dim_ineq = C_ineq.n_rows;
  
  if(dim_Q < 0)
    dim_Q = Q.n_rows;
  
  // -------------------------------------------------------
  // sqp::misc::section("Generate slack variables", debug);
  // sqp::misc::section("and calculate initial values", debug);
  // useful for rather bad initial guesses
  // but slow in convergence
  // -------------------------------------------------------
  
  arma::vec slack_eq_p,slack_eq_n,slack_ineq;
  arma::uvec equalities_above, equalities_below,inequalities_invalid;
  {
    const arma::vec equalities = C_eq * x - t_eq;
    const arma::vec inequalities = C_ineq * x - t_ineq;
    
    equalities_above = arma::find( equalities>tol);
    equalities_below = arma::find(-equalities>tol);
    inequalities_invalid = arma::find(inequalities>tol);
    
    if(all_slack)
    {
      slack_ineq.set_size(dim_ineq);
      
      slack_ineq.zeros();
      
      slack_eq_p.set_size(dim_eq);
      slack_eq_n.set_size(dim_eq);
      slack_eq_p.zeros();
      slack_eq_n.zeros();
      
      slack_eq_n.rows(equalities_above) = equalities.rows(equalities_above);
      slack_eq_p.rows(equalities_below) = -equalities.rows(equalities_below);
      
      slack_ineq.rows(inequalities_invalid) = inequalities.rows(inequalities_invalid);
      
      equalities_above = arma::linspace<arma::uvec>(0,dim_eq-1,dim_eq);
      equalities_below = arma::linspace<arma::uvec>(0,dim_eq-1,dim_eq);
      inequalities_invalid = arma::linspace<arma::uvec>(0,dim_ineq-1,dim_ineq);
    } else
    {
      if(equalities_above.n_elem>0)
        slack_eq_n = equalities.rows(equalities_above);
      if(equalities_below.n_elem>0)
        slack_eq_p = -equalities.rows(equalities_below);      
      if(inequalities_invalid.n_elem>0)
        slack_ineq = inequalities.rows(inequalities_invalid);
    }
  }
  
  // -------------------------------------------------------
  // sqp::misc::section("Augment sqp matrices", debug);
  // for slack variables
  // -------------------------------------------------------
  const unsigned dim_aug =  slack_eq_p.n_rows+slack_eq_n.n_rows+slack_ineq.n_rows;
  
  C_eq.resize(dim_eq,
              dim_Q+dim_aug); 
  
  // Equality slack variable addition to equality constraints
  if(equalities_below.n_elem>0)
  {
    const arma::uvec col_pos = arma::linspace<arma::uvec>(0,equalities_below.n_rows-1,equalities_below.n_rows);
    C_eq(equalities_below,col_pos+dim_Q) = arma::eye(equalities_below.n_rows,equalities_below.n_rows);
  }
  
  // Equality slack variable substraction from equality constraints
  if(equalities_above.n_elem>0)
  {    
    const arma::uvec col_pos = arma::linspace<arma::uvec>(0,equalities_above.n_rows-1,equalities_above.n_rows);
    C_eq(equalities_above,col_pos+slack_eq_p.n_rows+dim_Q) = -arma::eye(equalities_above.n_rows,equalities_above.n_rows);
  }
  
  // Rcpp::Rcout << "Adding " << dim_aug << " variables! " << "dim_ineq is " << dim_ineq << "\n";
  // Rcpp::Rcout << "C_ineq: " << C_ineq.n_rows << " x " << C_ineq.n_cols << "\n";
  
  C_ineq.resize(dim_ineq+dim_aug,
                dim_Q+dim_aug); 
  
  
  // Rcpp::Rcout << "C_ineq: " << C_ineq.n_rows << " x " << C_ineq.n_cols << "\n";
  
  // Slack variables must be >= 0
  if(dim_aug>0)
    C_ineq(dim_ineq,dim_Q,arma::size(dim_aug,dim_aug)).diag().fill(-1);
  // Substract inequality slack variables from inequality constraints
  
  if(inequalities_invalid.n_elem>0)
  {
    const arma::uvec col_pos = arma::linspace<arma::uvec>(0,inequalities_invalid.n_rows-1,inequalities_invalid.n_rows);
    
    C_ineq(inequalities_invalid,
           col_pos+slack_eq_p.n_rows+slack_eq_n.n_rows+dim_Q) = -arma::eye(inequalities_invalid.n_rows,inequalities_invalid.n_rows);
    
  }
  
  t_ineq.resize(dim_ineq + dim_aug);
  
  Q.resize(dim_Q+dim_aug,dim_Q+dim_aug);
  
  // Add penalty for slack variables
  l.resize(l.n_rows+dim_aug);
  if(dim_aug>0)
    l.subvec(dim_Q,l.n_rows-1).fill(penalty);
  
  
  x = arma::join_cols(arma::join_cols(x,slack_eq_p),arma::join_cols(slack_eq_n,slack_ineq));
  // -------------------------------------------------------
  // sqp::misc::section("Solve sqp", debug);
  // for slack variables
  // -------------------------------------------------------  
  
  sqp::solvers::result output = sqp::solvers::inequalities(x,
                                                         Q,
                                                         C_eq,
                                                         C_ineq,
                                                         l,
                                                         t_eq,
                                                         t_ineq,
                                                         tol,
                                                         max_iter,
                                                         C_eq.n_rows,
                                                         C_ineq.n_rows,
                                                         Q.n_rows,
                                                         solver,
                                                         fast,
                                                         debug);
  
  output.slack_eq_positive.set_size(dim_eq);
  output.slack_eq_negative.set_size(dim_eq);
  output.slack_ineq.set_size(dim_ineq);
  //
  output.lagrange_slack_eq_positive.set_size(dim_eq);
  output.lagrange_slack_eq_negative.set_size(dim_eq);
  output.lagrange_slack_ineq.set_size(dim_ineq);
  //
  output.slack_eq_positive.zeros();
  output.slack_eq_negative.zeros();
  output.slack_ineq.zeros();
  //
  output.lagrange_slack_eq_positive.zeros();
  output.lagrange_slack_eq_negative.zeros();
  output.lagrange_slack_ineq.zeros();
  
  if(equalities_below.n_rows>0)
  {
    output.slack_eq_positive.rows(equalities_below) = output.x.subvec(dim_Q,dim_Q+slack_eq_p.n_rows-1);
    output.lagrange_slack_eq_positive.rows(equalities_below) = output.lagrange_ineq.subvec(dim_ineq,dim_ineq+slack_eq_p.n_rows-1);
  }
  
  if(equalities_above.n_rows>0)
  {
    output.slack_eq_negative.rows(equalities_above) = output.x.subvec(dim_Q+slack_eq_p.n_rows,dim_Q+slack_eq_p.n_rows+slack_eq_n.n_rows-1);
    output.lagrange_slack_eq_negative.rows(equalities_above) = output.lagrange_ineq.subvec(dim_ineq+slack_eq_p.n_rows,dim_ineq+slack_eq_p.n_rows+slack_eq_n.n_rows-1);
  }
  
  if(inequalities_invalid.n_rows>0)
  {
    output.slack_ineq.rows(inequalities_invalid) = output.x.subvec(dim_Q+slack_eq_p.n_rows+slack_eq_n.n_rows,dim_Q+dim_aug-1);
    output.lagrange_slack_ineq.rows(inequalities_invalid) = output.lagrange_ineq.subvec(dim_ineq+slack_eq_p.n_rows+slack_eq_n.n_rows, dim_ineq+dim_aug-1);
  }
  
  if(dim_ineq>0)
  {
    output.lagrange_ineq = output.lagrange_ineq.subvec(0,dim_ineq-1);
  } else
  {
    output.lagrange_ineq.zeros(0);
  }
  
  output.x = output.x.subvec(0,dim_Q-1);
  
  return(output);
}

inline sqp::solvers::result slacked(
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
  // -------------------------------------------------------
  sqp::misc::section("Slacked sqp (sparse)", debug);
  // -------------------------------------------------------
  
  if(dim_eq < 0)
    dim_eq = C_eq.n_rows;
  
  if(dim_ineq < 0)
    dim_ineq = C_ineq.n_rows;
  
  if(dim_Q < 0)
    dim_Q = Q.n_rows;
  
  // -------------------------------------------------------
  sqp::misc::section("Check dimensions", debug);
  // -------------------------------------------------------
  
  if((unsigned)dim_Q != l.n_rows)
    Rcpp::stop(std::string("sqp error:\nDimensions of quadratic multiplier (")
                 + std::to_string(dim_Q) 
                 + std::string(") do not match with that of linear multiplier (")
                 + std::to_string(l.n_rows) 
                 + std::string(")")
    );
  
  if((unsigned)dim_Q != C_eq.n_cols)
    Rcpp::stop(std::string("sqp error:\nDimensions of quadratic multiplier (")
                 + std::to_string(dim_Q) 
                 + std::string(") do not match with that of equality constraint multiplier (")
                 + std::to_string(C_eq.n_cols) 
                 + std::string(")")
    );
  
  if((unsigned)dim_Q != C_ineq.n_cols)
    Rcpp::stop(std::string("sqp error:\nDimensions of quadratic multiplier (")
                 + std::to_string(dim_Q) 
                 + std::string(") do not match with that of inequality constraint multiplier (")
                 + std::to_string(C_ineq.n_cols) 
                 + std::string(")")
    );
  
  // -------------------------------------------------------
  // sqp::misc::section("Generate slack variables", debug);
  // sqp::misc::section("and calculate initial values", debug);
  // useful for rather bad initial guesses
  // but slow in convergence
  // -------------------------------------------------------
  
  arma::vec slack_eq_p,slack_eq_n,slack_ineq;
  arma::uvec equalities_above, equalities_below,inequalities_invalid;
  {
    // -------------------------------------------------------
    sqp::misc::section("Check constraints", debug);
    // -------------------------------------------------------
    const arma::vec equalities = C_eq * x - t_eq;
    const arma::vec inequalities = C_ineq * x - t_ineq;
    
    equalities_above = arma::find( equalities>tol);
    equalities_below = arma::find(-equalities>tol);
    inequalities_invalid = arma::find(inequalities>tol);
    
    // Rcpp::Rcout << "equalities above:\n " << equalities.rows(equalities_above).t() << "\n";
    // Rcpp::Rcout << "equalities below:\n " << equalities.rows(equalities_below).t() << "\n";
    // Rcpp::Rcout << "inequalities:\n " << inequalities.rows(inequalities_invalid).t() << "\n";
    
    if(all_slack)
    {
      slack_ineq.zeros(dim_ineq);
      slack_eq_p.zeros(dim_eq);
      slack_eq_n.zeros(dim_eq);
      
      slack_eq_n.rows(equalities_above) = equalities.rows(equalities_above);
      slack_eq_p.rows(equalities_below) = -equalities.rows(equalities_below);
      
      slack_ineq.rows(inequalities_invalid) = inequalities.rows(inequalities_invalid);
      
      equalities_above = arma::linspace<arma::uvec>(0,dim_eq-1,dim_eq);
      equalities_below = arma::linspace<arma::uvec>(0,dim_eq-1,dim_eq);
      inequalities_invalid = arma::linspace<arma::uvec>(0,dim_ineq-1,dim_ineq);
      
    } else
    {
      if(equalities_above.n_elem>0)
        slack_eq_n = equalities.rows(equalities_above);
      if(equalities_below.n_elem>0)
        slack_eq_p = -equalities.rows(equalities_below);      
      if(inequalities_invalid.n_elem>0)
        slack_ineq = inequalities.rows(inequalities_invalid);
    }
  }
  // -------------------------------------------------------
  // sqp::misc::section("Augment sqp matrices", debug);
  // for slack variables
  // -------------------------------------------------------
  const unsigned dim_aug =  slack_eq_p.n_rows+slack_eq_n.n_rows+slack_ineq.n_rows;
  
  if(dim_aug>0)
    C_eq.resize(dim_eq,
                dim_Q+dim_aug);
  
  // Equality slack variable addition to equality constraints
  if(equalities_below.n_elem>0)
  {
    for (unsigned j = 0; j < equalities_below.n_rows; ++j){
      unsigned i = arma::as_scalar(equalities_below.row(j));
      C_eq(i,j+dim_Q) = 1;
    }
  }
  
  // Equality slack variable substraction from equality constraints
  if(equalities_above.n_elem>0)
  {
    for (unsigned j = 0; j < equalities_above.n_rows; ++j){
      unsigned i = arma::as_scalar(equalities_above.row(j));
      C_eq(i,j+slack_eq_p.n_rows+dim_Q) = -1;
    }
    
  }
  
  C_ineq.resize(dim_ineq+dim_aug,
                dim_Q+dim_aug);
  
  // Slack variables must be >= 0
  if(dim_aug>0)
    C_ineq(dim_ineq,dim_Q,arma::size(dim_aug,dim_aug)) = -arma::speye<arma::sp_mat>(dim_aug,dim_aug);//
  
  // Inequality slack variable substraction from inequality constraints
  if(inequalities_invalid.n_elem>0)
  {
    for (unsigned j = 0; j < inequalities_invalid.n_rows; ++j){
      unsigned i = arma::as_scalar(inequalities_invalid.row(j));
      C_ineq(i,j+slack_eq_p.n_rows+slack_eq_n.n_rows+dim_Q) = -1;//arma::eye(equalities_below.n_rows,equalities_below.n_rows);
    }
  }
  t_ineq.resize(dim_ineq + dim_aug);
  Q.resize(dim_Q+dim_aug,dim_Q+dim_aug);
  
  
  // Add penalty for slack variables
  l.resize(l.n_rows+dim_aug);
  if(dim_aug>0)
  {
    l.subvec(dim_Q,l.n_rows-1).fill(penalty);
    t_ineq.subvec(dim_ineq,t_ineq.n_rows-1).zeros();
  }
  
  x = arma::join_cols(arma::join_cols(x,slack_eq_p),arma::join_cols(slack_eq_n,slack_ineq));
  // -------------------------------------------------------
  // sqp::misc::section("Solve sqp", debug);
  // for slack variables
  // -------------------------------------------------------
  
  // arma::mat Qx(Q);
  // arma::mat C_eqx(C_eq);
  // arma::mat C_ineqx(C_ineq);
  // // 
  // const arma::uvec index_Q = arma::find(arma::max(Qx,1) == 0  && arma::min(Qx,1) == 0);
  // const arma::uvec index_C_eq = arma::find(arma::max(C_eqx,1) == 0  && arma::min(C_eqx,1) == 0);
  // const arma::uvec index_C_ineq = arma::find(arma::max(C_ineqx,1) == 0  && arma::min(C_ineqx,1) == 0);
  // 
  // const arma::uvec unused_Q = arma::find(arma::max(Qx,0) == 0  && arma::min(Qx,0) == 0);
  // const arma::uvec unused_eq = arma::find(arma::max(C_eqx,0) == 0  && arma::min(C_eqx,0) == 0);
  // const arma::uvec unused_ineq = arma::find(arma::max(C_ineqx,0) == 0  && arma::min(C_ineqx,0) == 0);
  // 
  // const arma::uvec unused = arma::intersect(unused_Q,arma::intersect(unused_eq,unused_ineq));
  // 
  // Rcpp::Rcout << "index_Q: " << index_Q.t() << "\n";
  // Rcpp::Rcout << "index_C_eq: " << index_C_eq.t() << "\n";
  // Rcpp::Rcout << "index_C_ineq: " << index_C_ineq.t() << "\n";
  // 
  // Rcpp::Rcout << "Variables not used for distance: " << unused_Q.t() << "\n";
  // Rcpp::Rcout << "Variables not used for equalities: " << unused_eq.t() << "\n";
  // Rcpp::Rcout << "Variables not used for inequalities: " << unused_ineq.t() << "\n";
  // Rcpp::Rcout << "Variables not used: " << unused.t() << "\n";
  // 
  // Rcpp::Rcout << "-----------------------------------------------\n";
  
  sqp::solvers::result output = sqp::solvers::inequalities(x,Q,C_eq,C_ineq,l,t_eq,t_ineq,tol,max_iter,C_eq.n_rows,C_ineq.n_rows,Q.n_rows,debug);
  
  output.slack_eq_positive.set_size(dim_eq);
  output.slack_eq_negative.set_size(dim_eq);
  output.slack_ineq.set_size(dim_ineq);
  //
  output.lagrange_slack_eq_positive.set_size(dim_eq);
  output.lagrange_slack_eq_negative.set_size(dim_eq);
  output.lagrange_slack_ineq.set_size(dim_ineq);
  //
  output.slack_eq_positive.zeros();
  output.slack_eq_negative.zeros();
  output.slack_ineq.zeros();
  //
  output.lagrange_slack_eq_positive.zeros();
  output.lagrange_slack_eq_negative.zeros();
  output.lagrange_slack_ineq.zeros();
  
  if(equalities_below.n_rows>0)
  {
    output.slack_eq_positive.rows(equalities_below) = output.x.subvec(dim_Q,dim_Q+slack_eq_p.n_rows-1);
    output.lagrange_slack_eq_positive.rows(equalities_below) = output.lagrange_ineq.subvec(dim_ineq,dim_ineq+slack_eq_p.n_rows-1);
  }
  
  if(equalities_above.n_rows>0)
  {
    output.slack_eq_negative.rows(equalities_above) = output.x.subvec(dim_Q+slack_eq_p.n_rows,dim_Q+slack_eq_p.n_rows+slack_eq_n.n_rows-1);
    output.lagrange_slack_eq_negative.rows(equalities_above) = output.lagrange_ineq.subvec(dim_ineq+slack_eq_p.n_rows,dim_ineq+slack_eq_p.n_rows+slack_eq_n.n_rows-1);
  }
  
  if(inequalities_invalid.n_rows>0)
  {
    output.slack_ineq.rows(inequalities_invalid) = 
      output.x.subvec(dim_Q+slack_eq_p.n_rows+slack_eq_n.n_rows,dim_Q+dim_aug-1);
    
    output.lagrange_slack_ineq.rows(inequalities_invalid) = 
      output.lagrange_ineq.subvec(dim_ineq+slack_eq_p.n_rows+slack_eq_n.n_rows, dim_ineq+dim_aug-1);
  }
  
  if(dim_ineq>0)
  {
    output.lagrange_ineq = output.lagrange_ineq.subvec(0,dim_ineq-1);
  }
  else
  {
    output.lagrange_ineq.zeros(0);
  }
  output.x = output.x.subvec(0,dim_Q-1);
  
  return(output);
}



}
}

#endif                                                 // end of include guard for 'solvers/quadratic_slack.h'  
