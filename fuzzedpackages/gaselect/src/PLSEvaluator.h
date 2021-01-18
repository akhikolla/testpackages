//
//  PLSEvaluator.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 03.05.2013.
//
//

#ifndef GenAlgPLS_PLSEvaluator_h
#define GenAlgPLS_PLSEvaluator_h

#include "config.h"

#include <memory>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <RcppArmadillo.h>

#include "RNG.h"
#include "Evaluator.h"
#include "Chromosome.h"
#include "PLS.h"

class PLSEvaluator : public Evaluator {
public:
	PLSEvaluator(std::unique_ptr<PLS> pls, uint16_t numReplications, uint16_t maxNComp, const std::vector<uint32_t> &seed,
               VerbosityLevel verbosity, uint16_t innerSegments, uint16_t outerSegments = 1, double testSetSize = 0.0,
               double sdfact = 1.0);

	~PLSEvaluator() {}

	double evaluate(Chromosome &ch) {
		arma::uvec columnSubset = ch.toColumnSubset();
		double fitness = this->evaluate(columnSubset);
		ch.setFitness(fitness);
		return fitness;
	};

	double evaluate(arma::uvec &columnSubset);

	std::vector<arma::uvec> getSegmentation() const {
		return this->segmentation;
	}

	Evaluator* clone() const;

#ifdef ENABLE_DEBUG_VERBOSITY
	static uint32_t counter;
#endif

private:
	const uint16_t numReplications;
	const uint16_t outerSegments;
	const uint16_t innerSegments;
	const double sdfact;
	const arma::uword nrows;

	std::unique_ptr<PLS> pls;
	uint16_t maxNComp;
	std::vector<arma::uvec> segmentation;

	PLSEvaluator(const PLSEvaluator &other);

	/**
	 * Estimate the SEP
	 */
	double estSEP(uint16_t maxNComp);

	void initSegmentation(double testSetSize, const std::vector<uint32_t> &seed);

};


#endif
