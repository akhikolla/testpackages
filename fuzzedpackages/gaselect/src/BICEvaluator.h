//
//  BICEvaluator.h
//  gaselect
//
//  Created by David Kepplinger on 09.08.2014.
//
//

#ifndef gaselect_BICEvaluator_h
#define gaselect_BICEvaluator_h

#include "config.h"

#include <stdexcept>
#include <vector>
#include <memory>
#include <RcppArmadillo.h>

#include "RNG.h"
#include "Evaluator.h"
#include "Chromosome.h"
#include "PLS.h"

class BICEvaluator : public Evaluator {
public:
	enum Statistic {
		BIC = 0,
		AIC = 1,
		ADJ_R2 = 2,
		R2 = 3
	};

	BICEvaluator(std::unique_ptr<PLS> pls, uint16_t maxNComp, const std::vector<uint32_t> &seed, VerbosityLevel verbosity,
		uint16_t numSegments = 7, Statistic stat = BIC, double sdfact = 1.0);

	~BICEvaluator() {};

	double evaluate(Chromosome &ch) {
		arma::uvec columnSubset = ch.toColumnSubset();
		double fitness = this->evaluate(columnSubset);
		ch.setFitness(fitness);
		return fitness;
	};

	std::vector<arma::uvec> getSegmentation() const {
		return this->segmentation;
	}

	double evaluate(arma::uvec &columnSubset);

	Evaluator* clone() const;

private:
	const uint16_t numSegments;
	const arma::uword nrows;
	const double sdfact;
	const BICEvaluator::Statistic stat;

	std::unique_ptr<PLS> pls;
	uint16_t maxNComp;
	std::vector<arma::uvec> segmentation;
	double r2denom;

	BICEvaluator(const BICEvaluator &other);

	/**
	 * Estimate the SEP
	 */
	double getRSS(uint16_t maxNComp);

	void initSegmentation(const std::vector<uint32_t> &seed);
};

#endif
