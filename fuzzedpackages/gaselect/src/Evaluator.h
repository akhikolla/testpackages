//
//  Evaluator.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 15.04.2013.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef GenAlgPLS_Evaluator_h
#define GenAlgPLS_Evaluator_h

#include "config.h"
#include <vector>
#include <RcppArmadillo.h>
#include "Control.h"
#include "Chromosome.h"
#include "RNG.h"

class Evaluator {
public:
	class EvaluatorException : public std::runtime_error {
	public:
		EvaluatorException(const char* what) : std::runtime_error(what) {};
		virtual ~EvaluatorException() throw() {};
	};

	Evaluator(const VerbosityLevel verbosity) : verbosity(verbosity) {}
	virtual ~Evaluator() {};

	virtual double evaluate(arma::uvec &columnSubset) = 0;
	virtual double evaluate(Chromosome &ch) = 0;
	
	virtual Evaluator* clone() const = 0;

	virtual std::vector<arma::uvec> getSegmentation() const {
		return std::vector<arma::uvec>();
	}
protected:
	const VerbosityLevel verbosity;
};


#endif
