#include "optimization_problem.h"

// [[Rcpp::export]]
SEXP rcpp_new_optimization_problem(std::size_t nrow = 1000000,
                                   std::size_t ncol = 1000000,
                                   std::size_t ncell = 100000) {
  OPTIMIZATIONPROBLEM* x = new OPTIMIZATIONPROBLEM(nrow, ncol, ncell);
  Rcpp::XPtr<OPTIMIZATIONPROBLEM> ptr =
    Rcpp::XPtr<OPTIMIZATIONPROBLEM>(x, true);
  return(ptr);
}

// [[Rcpp::export]]
SEXP rcpp_predefined_optimization_problem(Rcpp::List l) {
  std::string modelsense = Rcpp::as<std::string>(l["modelsense"]);
  std::size_t number_of_projects =
    Rcpp::as<std::size_t>(l["number_of_projects"]);
  std::size_t number_of_actions =
    Rcpp::as<std::size_t>(l["number_of_actions"]);
  std::size_t number_of_features =
    Rcpp::as<std::size_t>(l["number_of_features"]);
  std::size_t number_of_branches =
    Rcpp::as<std::size_t>(l["number_of_branches"]);
  std::vector<std::size_t> A_i = Rcpp::as<std::vector<std::size_t>>(l["A_i"]);
  std::vector<std::size_t> A_j = Rcpp::as<std::vector<std::size_t>>(l["A_j"]);
  std::vector<double> A_x = Rcpp::as<std::vector<double>>(l["A_x"]);
  std::vector<double> obj = Rcpp::as<std::vector<double>>(l["obj"]);
  Rcpp::List pwlobj = Rcpp::as<Rcpp::List>(l["pwlobj"]);
  std::vector<double> lb = Rcpp::as<std::vector<double>>(l["lb"]);
  std::vector<double> ub = Rcpp::as<std::vector<double>>(l["ub"]);
  std::vector<double> rhs = Rcpp::as<std::vector<double>>(l["rhs"]);
  std::vector<std::string> sense =
    Rcpp::as<std::vector<std::string>>(l["sense"]);
  std::vector<std::string> vtype =
    Rcpp::as<std::vector<std::string>>(l["vtype"]);
  std::vector<std::string> row_ids =
    Rcpp::as<std::vector<std::string>>(l["row_ids"]);
  std::vector<std::string> col_ids =
    Rcpp::as<std::vector<std::string>>(l["col_ids"]);
  OPTIMIZATIONPROBLEM* x = new OPTIMIZATIONPROBLEM(modelsense,
   number_of_projects, number_of_actions, number_of_features,
   number_of_branches, A_i, A_j, A_x, obj, pwlobj, lb, ub, rhs, sense, vtype,
   row_ids, col_ids);
  Rcpp::XPtr<OPTIMIZATIONPROBLEM> ptr =
    Rcpp::XPtr<OPTIMIZATIONPROBLEM>(x, true);
  return(ptr);
}

// [[Rcpp::export]]
Rcpp::List rcpp_optimization_problem_as_list(SEXP x) {
  // initialization
  Rcpp::XPtr<OPTIMIZATIONPROBLEM> ptr =
    Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x);
  // create list
  return Rcpp::List::create(
    Rcpp::Named("modelsense") = ptr->_modelsense,
    Rcpp::Named("number_of_projects") = ptr->_number_of_projects,
    Rcpp::Named("number_of_actions") = ptr->_number_of_actions,
    Rcpp::Named("number_of_features") = ptr->_number_of_features,
    Rcpp::Named("number_of_branches") = ptr->_number_of_branches,
    Rcpp::Named("A_i") = Rcpp::IntegerVector(ptr->_A_i.begin(),
                                             ptr->_A_i.end()),
    Rcpp::Named("A_j") = Rcpp::IntegerVector(ptr->_A_j.begin(),
                                             ptr->_A_j.end()),
    Rcpp::Named("A_x") = ptr->_A_x,
    Rcpp::Named("obj") = ptr->_obj,
    Rcpp::Named("pwlobj") = ptr->_pwlobj,
    Rcpp::Named("lb") = ptr->_lb,
    Rcpp::Named("ub") = ptr->_ub,
    Rcpp::Named("rhs") = ptr->_rhs,
    Rcpp::Named("sense") = ptr->_sense,
    Rcpp::Named("vtype") = ptr->_vtype,
    Rcpp::Named("row_ids") = ptr->_row_ids,
    Rcpp::Named("col_ids") = ptr->_col_ids);
}

// [[Rcpp::export]]
std::size_t rcpp_get_optimization_problem_ncol(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->ncol());
}

// [[Rcpp::export]]
std::size_t rcpp_get_optimization_problem_nrow(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->nrow());
}

// [[Rcpp::export]]
std::size_t rcpp_get_optimization_problem_ncell(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->ncell());
}

// [[Rcpp::export]]
Rcpp::List rcpp_get_optimization_problem_A(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->A());
}

// [[Rcpp::export]]
std::string rcpp_get_optimization_problem_modelsense(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_modelsense);
}

// [[Rcpp::export]]
std::size_t rcpp_get_optimization_problem_number_of_projects(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_number_of_projects);
}

// [[Rcpp::export]]
std::size_t rcpp_get_optimization_problem_number_of_actions(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_number_of_actions);
}

// [[Rcpp::export]]
std::size_t rcpp_get_optimization_problem_number_of_features(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_number_of_features);
}

// [[Rcpp::export]]
std::size_t rcpp_get_optimization_problem_number_of_branches(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_number_of_branches);
}

// [[Rcpp::export]]
std::vector<std::string> rcpp_get_optimization_problem_vtype(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_vtype);
}

// [[Rcpp::export]]
std::vector<double> rcpp_get_optimization_problem_obj(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_obj);
}

// [[Rcpp::export]]
Rcpp::List rcpp_get_optimization_problem_pwlobj(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_pwlobj);
}

// [[Rcpp::export]]
std::vector<double> rcpp_get_optimization_problem_rhs(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_rhs);
}

// [[Rcpp::export]]
std::vector<std::string> rcpp_get_optimization_problem_sense(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_sense);
}

// [[Rcpp::export]]
std::vector<double> rcpp_get_optimization_problem_lb(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_lb);
}

// [[Rcpp::export]]
std::vector<double> rcpp_get_optimization_problem_ub(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_ub);
}

// [[Rcpp::export]]
std::vector<std::string> rcpp_get_optimization_problem_col_ids(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_col_ids);
}

// [[Rcpp::export]]
std::vector<std::string> rcpp_get_optimization_problem_row_ids(SEXP x) {
  return(Rcpp::as<Rcpp::XPtr<OPTIMIZATIONPROBLEM>>(x)->_row_ids);
}
