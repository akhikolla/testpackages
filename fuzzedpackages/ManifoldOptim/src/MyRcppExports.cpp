// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "ManifoldOptimProblem.h"

using namespace Rcpp;

// The Problem object is passed from R, but it actually represents a C++ class.
// We want to take the Problem object back to C++ now. The solution below was
// suggested in a thread:
// See: http://rcpp-devel.r-forge.r-project.narkive.com/5MzKnjRd/
//	rcpp-modules-passing-object-pointer-back-to-c

Rcpp::List ManifoldOptim(const arma::vec& initX, const arma::mat& initH,
	ManifoldOptimProblem& prob, const Rcpp::List& maniDefn, const Rcpp::List& maniParams,
	const Rcpp::List& solverParams, const Rcpp::List& derivParams,
	const std::string&  method, bool hasHHR);
RcppExport SEXP manifold_optim(SEXP sexp_initX, SEXP sexp_initH,
	SEXP sexp_prob, SEXP sexp_maniDefn, SEXP sexp_maniParams,
	SEXP sexp_solverParams, SEXP sexp_derivParams,
	SEXP sexp_method, SEXP sexp_hasHHR)
{
BEGIN_RCPP
	Rcpp::RObject __result;
	Rcpp::RNGScope __rngScope;
	Rcpp::traits::input_parameter<const arma::vec&>::type initX(sexp_initX);
	Rcpp::traits::input_parameter<const arma::mat&>::type initH(sexp_initH);

	std::string rtypename("ManifoldOptimProblem");
	Rcpp::S4 s4obj(sexp_prob);
	Rcpp::Environment env(s4obj);
	Rcpp::XPtr<ManifoldOptimProblem> xptr( env.get(".pointer") );
	ManifoldOptimProblem* prob_ptr = static_cast<ManifoldOptimProblem*> (R_ExternalPtrAddr(xptr));

	Rcpp::traits::input_parameter<const List&>::type maniDefn(sexp_maniDefn);
	Rcpp::traits::input_parameter<const List&>::type maniParams(sexp_maniParams);
	Rcpp::traits::input_parameter<const List&>::type solverParams(sexp_solverParams);
	Rcpp::traits::input_parameter<const List&>::type derivParams(sexp_derivParams);
	Rcpp::traits::input_parameter<std::string>::type method(sexp_method);
	Rcpp::traits::input_parameter<bool>::type hasHHR(sexp_hasHHR);
	__result = Rcpp::wrap(ManifoldOptim(initX, initH, *prob_ptr, maniDefn, maniParams, solverParams, derivParams, method, hasHHR));
	return __result;
END_RCPP
}

