#include "RProblem.h"

/*
 * The ls Rcpp::Function is used as a default value, indicating that
 * no gradient or Hessian has been set by the user. There's no special reason
 * we used ls here, only tha it would never be used to compute a gradient or Hessian.
 */
RProblem::RProblem(const Rcpp::Function& objFun)
: ManifoldOptimProblem(), m_objFun(objFun), 
  m_gradFun(Rcpp::Function("ls")),
  m_hessEtaFun(Rcpp::Function("ls")),
  m_defaultFun(Rcpp::Function("ls"))
{
}

RProblem::RProblem(const Rcpp::Function& objFun, const Rcpp::Function& gradFun)
: ManifoldOptimProblem(), m_objFun(objFun), m_gradFun(gradFun),
  m_hessEtaFun(Rcpp::Function("ls")),
  m_defaultFun(Rcpp::Function("ls"))
{
}

RProblem::RProblem(const Rcpp::Function& objFun, const Rcpp::Function& gradFun,
	const Rcpp::Function& hessEtaFun)
: ManifoldOptimProblem(), m_objFun(objFun), m_gradFun(gradFun),
  m_hessEtaFun(hessEtaFun), m_defaultFun(Rcpp::Function("ls"))
{
}

double RProblem::objFun(const arma::vec &X) const
{
	return Rcpp::as<double>(m_objFun(X));
}

arma::vec RProblem::gradFun(const arma::vec &X) const
{
	// If no Rcpp::Function was specified, call the base class gradFun
	if (m_gradFun == m_defaultFun) {
		return ManifoldOptimProblem::gradFun(X);
	}
	return Rcpp::as<arma::vec>(m_gradFun(X));
}

arma::vec RProblem::hessEtaFun(const arma::vec &X, const arma::vec &eta) const
{
	// If no Rcpp::Function was specified, call the base class hessEtaFun
	if (m_hessEtaFun == m_defaultFun) {
		return ManifoldOptimProblem::hessEtaFun(X, eta);
	}
	return Rcpp::as<arma::vec>(m_hessEtaFun(X, eta));
}

