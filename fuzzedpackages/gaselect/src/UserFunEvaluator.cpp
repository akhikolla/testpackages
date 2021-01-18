//
//  UserFunEvaluator.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 03.05.2013.
//
//

#include "config.h"
#include <RcppArmadillo.h>
#include "UserFunEvaluator.h"

double UserFunEvaluator::evaluate(Chromosome &ch) {
	SEXP rawFitness = this->userFun(Rcpp::wrap(ch.toLogicalVector()));
	if(!Rf_isNumeric(rawFitness)) {
		throw Rcpp::exception("Evaluation function has to return a numeric value", __FILE__, __LINE__);
	}

	double fitness = Rcpp::as<double>(rawFitness);
	ch.setFitness(fitness);
	return fitness;
}

double UserFunEvaluator::evaluate(arma::uvec &columnSubset) {
	Rcpp::LogicalVector logVec(columnSubset.n_elem, false);

	for(arma::uword i = 0; i < columnSubset.n_elem; ++i) {
		logVec[columnSubset[i]] = true;
	}
	
	SEXP rawFitness = this->userFun(Rcpp::wrap(logVec));
	if(!Rf_isNumeric(rawFitness)) {
		throw Evaluator::EvaluatorException("The evaluation function did not return a numeric value.");
	}
	
	double fitness = Rcpp::as<double>(rawFitness);
	return fitness;
}
