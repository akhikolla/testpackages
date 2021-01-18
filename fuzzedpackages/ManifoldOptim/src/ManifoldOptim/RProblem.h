#ifndef R_PROBLEM_H
#define R_PROBLEM_H

#include <RcppArmadillo.h>
#include "ManifoldOptimProblem.h"
#include <cstring>
#include <def.h>
#include <assert.h>
#include <iostream>

using namespace Rcpp;

class RProblem : public ManifoldOptimProblem
{

public:
	RProblem(const Rcpp::Function& objFun);
	RProblem(const Rcpp::Function& objFun, const Rcpp::Function& gradFun);
	RProblem(const Rcpp::Function& objFun, const Rcpp::Function& gradFun,
		const Rcpp::Function& hessEtaFun);

	double objFun(const arma::vec& X) const;
	arma::vec gradFun(const arma::vec& X) const;
	arma::vec hessEtaFun(const arma::vec& X, const arma::vec& eta) const;

private:
	Rcpp::Function m_objFun;
	Rcpp::Function m_gradFun;
	Rcpp::Function m_hessEtaFun;
	Rcpp::Function m_defaultFun;
};

#endif

