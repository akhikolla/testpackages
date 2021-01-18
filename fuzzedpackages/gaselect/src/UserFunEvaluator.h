//
//  UserFunEvaluator.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 03.05.2013.
//
//

#ifndef GenAlgPLS_UserFunEvaluator_h
#define GenAlgPLS_UserFunEvaluator_h

#include "config.h"

#include <stdexcept>
#include <RcppArmadillo.h>
#include "Evaluator.h"
#include "Chromosome.h"

class UserFunEvaluator : public Evaluator {
public:
	UserFunEvaluator(Rcpp::Function const &userFun, const VerbosityLevel &verbosity) : Evaluator(verbosity), userFun(userFun) {};
//	~UserFunEvaluator();

	double evaluate(Chromosome &ch);
	double evaluate(arma::uvec &columnSubset);
	/**
	 * The UserFunEvaluator can not be cloned!!
	 * It throws an std::logic_error if called
	 */
	Evaluator* clone() const { throw std::logic_error("A user specified evaluation function can not be cloned!"); }
private:
	const Rcpp::Function userFun;
};

#endif
