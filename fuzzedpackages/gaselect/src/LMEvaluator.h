//
//  LMEvaluator.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 28.08.2013.
//
//

#ifndef GenAlgPLS_LMEvaluator_h
#define GenAlgPLS_LMEvaluator_h

#include "config.h"

#include <exception>
#include <RcppArmadillo.h>
#include "Evaluator.h"
#include "Chromosome.h"

class LMEvaluator : public Evaluator {
public:
	enum Statistic {
		BIC = 0,
		AIC = 1,
		ADJ_R2 = 2,
		R2 = 3
	};
	
	LMEvaluator(const arma::mat &X, const arma::colvec &y, const LMEvaluator::Statistic statistic, const VerbosityLevel &verbosity, const bool addIntercept = true);
	//	~LMEvaluator();
	
	double evaluate(arma::uvec &columnSubset);
	
	double evaluate(Chromosome &ch) {
		arma::uvec columnSubset = ch.toColumnSubset();
		double fitness = this->evaluate(columnSubset);
		ch.setFitness(fitness);
		return fitness;
	};
	/**
	 * The UserFunEvaluator can not be cloned!!
	 * It throws an std::logic_error if called
	 */
	Evaluator* clone() const { return new LMEvaluator(this->Xdesign, this->y, this->statistic, this->verbosity, false); }
private:
	const arma::colvec y;
	const LMEvaluator::Statistic statistic;

	arma::mat Xdesign; // X matrix with a 1 column in front
	double r2denom;
};

#endif
