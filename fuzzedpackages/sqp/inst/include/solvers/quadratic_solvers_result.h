#ifndef __sqp_solvers_quadratic_solvers_result_included__        // if include guard for 'solvers/quadratic_solvers_result.h' is undefined
#define __sqp_solvers_quadratic_solvers_result_included__        // define include guard for 'solvers/quadratic_solvers_result.h'  

#include "sqp.h"

namespace sqp {
namespace solvers{



class result  {
public:
  arma::vec x{};
  arma::vec lagrange_eq{};
  arma::vec lagrange_ineq{};
  arma::vec slack_eq_positive{};
  arma::vec slack_eq_negative{};
  arma::vec slack_ineq{};
  arma::vec lagrange_slack_eq_positive{};
  arma::vec lagrange_slack_eq_negative{};
  arma::vec lagrange_slack_ineq{};
  
  void zeros()
  {
    x.zeros();
  };
  
  
  operator Rcpp::List() const { 
    return Rcpp::List::create(Rcpp::Named("x") = x,
                              Rcpp::Named("lagrange_eq") = lagrange_eq,
                              Rcpp::Named("lagrange_ineq") = lagrange_ineq,
                              Rcpp::Named("slack_eq_positive") = slack_eq_positive,
                              Rcpp::Named("slack_eq_negative") = slack_eq_negative,
                              Rcpp::Named("slack_ineq") = slack_ineq,
                              Rcpp::Named("lagrange_slack_eq_positive") = lagrange_slack_eq_positive,
                              Rcpp::Named("lagrange_slack_eq_negative") = lagrange_slack_eq_negative,
                              Rcpp::Named("lagrange_slack_ineq") = lagrange_slack_ineq
    );
  }
};



}
}

#endif                                                          // end of include guard for 'solvers/quadratic_solvers_result.h'  
