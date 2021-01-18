//
//  BICEvaluator.cpp
//  gaselect
//
//  Created by David Kepplinger on 09.08.2014.
//
//
#include "config.h"

#include <algorithm>
#include <memory>
#include "Logger.h"
#include "BICEvaluator.h"
#include "OnlineStddev.h"

#ifdef ENABLE_DEBUG_VERBOSITY
#define IF_DEBUG(expr) if(this->verbosity == DEBUG_EVAL || this->verbosity == DEBUG_ALL) { expr; }
#else
#define IF_DEBUG(expr)
#endif

BICEvaluator::BICEvaluator(std::unique_ptr<PLS> _pls, uint16_t _maxNComp, const std::vector<uint32_t> &seed,
                           VerbosityLevel _verbosity, uint16_t _numSegments, BICEvaluator::Statistic _stat,
                           double _sdfact) :
		Evaluator(_verbosity), numSegments(_numSegments),	nrows(_pls->getNumberOfObservations()),
		sdfact(_sdfact / sqrt((double) _numSegments)), stat(_stat), pls(std::move(_pls)), maxNComp(_maxNComp)
{
	if(pls->getNumberOfResponseVariables() > 1) {
		throw std::invalid_argument("PLS evaluator only available for models with 1 response variable");
	}

	if(numSegments < 2) {
		throw std::invalid_argument("For CV at least 2 segments are needed");
	}

//	this->r2denom = arma::sum(arma::square(this->pls->getY() - arma::mean(this->pls->getY())));
	this->r2denom = this->pls->getNumberOfObservations() * arma::var(this->pls->getY(), 1); // N * Var(Y)

	if(this->maxNComp <= 1) {
		this->maxNComp = this->nrows - 1;
	}

	this->initSegmentation(seed);
}

BICEvaluator::BICEvaluator(const BICEvaluator &other) :
	Evaluator(other.verbosity), numSegments(other.numSegments), nrows(other.nrows),
	sdfact(other.sdfact), stat(other.stat),
	maxNComp(other.maxNComp), segmentation(other.segmentation), r2denom(other.r2denom)
{
	this->pls = other.pls->clone();
}

double BICEvaluator::evaluate(arma::uvec &columnSubset) {
	if(columnSubset.n_elem == 0) {
		GAerr << GAerr.lock() << "Can not evaluate empty variable subset" << GAerr.unlock();
		throw Evaluator::EvaluatorException("Can not evaluate empty variable subset");
	}

	this->pls->viewSelectColumns(columnSubset);

	double ret = 0.0;
	double RSS = this->getRSS(((this->maxNComp < columnSubset.n_elem) ? this->maxNComp : columnSubset.n_elem));


	IF_DEBUG(GAout << "RSS:" << std::endl << RSS << std::endl)

	switch(this->stat) {
		case BIC: {
			ret = -(this->nrows * log(RSS / this->nrows) + columnSubset.n_elem * log((double) this->nrows));
			break;
		}
		case AIC: {
			ret = -(this->nrows * log(RSS / this->nrows) + columnSubset.n_elem * 2.0);
			break;
		}
		case ADJ_R2: {
			double r2 = 1 - (RSS / this->r2denom);
			ret = 1 - (((this->nrows - 1) * (1 - r2)) / (this->nrows - columnSubset.n_elem - 1));
			break;
		}
		case R2: {
			ret = 1 - (RSS / this->r2denom);
			break;
		}
		default:
			break;
	}

	return ret;
}

/**
 * Initialize the row-segmentation for each replication and all segmentations
 */
inline void BICEvaluator::initSegmentation(const std::vector<uint32_t> &seed) {
	RNG rng(seed);
	arma::uvec shuffledRowNumbers = ShuffledSet(this->nrows).shuffleAll(rng);
	arma::uword segmentLength = this->nrows / this->numSegments;
	arma::uword segmentLengthRem = this->nrows % this->numSegments;
	arma::uword segLen = 0, n = 0;


	this->segmentation.reserve(2 * this->numSegments);

	if((this->nrows - segmentLength - 2) < this->maxNComp) {
		this->maxNComp = this->nrows - segmentLength - 2;
	}

	IF_DEBUG(GAout << "Initialize segments with a segment length of " << segmentLength << std::endl);

	for(uint16_t i = 0; i < this->numSegments; ++i) {
		segLen = segmentLength;
		if(i < segmentLengthRem) {
			++segLen;
		}

		arma::uvec inSegment(this->nrows - segLen);

		if(n > 0) {
			inSegment.rows(0, n - 1) = shuffledRowNumbers.rows(0, n - 1);
		}

		if(n < this->nrows - segLen) {
			inSegment.rows(n, inSegment.n_elem - 1) = shuffledRowNumbers.rows(n + segLen, this->nrows - 1);
		}

		std::sort(inSegment.begin(), inSegment.end());

		/* First push back training set */
		this->segmentation.push_back(inSegment);
		IF_DEBUG(this->segmentation.back().t().raw_print(GAout, "Training set:"));
		/* Then add test set */
		this->segmentation.push_back(arma::sort(shuffledRowNumbers.rows(n, n + segLen - 1)));
		IF_DEBUG(this->segmentation.back().t().raw_print(GAout, "Test set:"));

		n += segLen;
	}
}

/**
 * Estimate the SEP for the given row ordering
 */
double BICEvaluator::getRSS(uint16_t maxNComp) {
	double RSS = 0.0;
	OnlineStddev trainMSEP(maxNComp);

	double cutoff;
	arma::uword minNComp = 0, optNComp = 0;

	arma::mat leftOutX;
	arma::vec leftOutY;

	uint16_t seg = 0, comp;
	std::vector<arma::uvec>::const_iterator segmentIter = this->segmentation.begin();

	try {
		/*
		 * Fit PLS models to predict the values in each segment once
		 */
		while(seg++ < this->numSegments) {
			/* Segmentation iterator currently points to the `fit` segment */
			this->pls->viewSelectRows(*(segmentIter));
			this->pls->fit(maxNComp);

			/* Increment segmentation iterator to point to the `predict` segment */
			++segmentIter;

			leftOutY = this->pls->getY().rows(*segmentIter);
			leftOutX = this->pls->getXColumnView().rows(*segmentIter);

			for(comp = 0; comp < maxNComp; ++comp) {
				trainMSEP.update(arma::mean(arma::square(leftOutY - this->pls->predict(leftOutX, comp + 1))), comp);
			}

			/* Increment segmentation iterator to point to the next `fit` segment */
			++segmentIter;
		}

		/*
		 * Find best number of components based on the RSS plus one standard deviation
		 */
		IF_DEBUG(
			GAout << "EVALUATOR: MSE and SD" << std::endl;
			for(uint16_t j = 0; j < maxNComp; ++j) {
				GAout << "\t (" << j + 1 << " comps.): " << trainMSEP.mean(j) << " +- " << trainMSEP.stddev(j) << std::endl;
			}
		)

		cutoff = trainMSEP.mean(0);
		minNComp = 0;

		for(comp = 1; comp < maxNComp; ++comp) {
			if(trainMSEP.mean(comp) < cutoff) {
				minNComp = comp;
				cutoff = trainMSEP.mean(comp);
			}
		}

		IF_DEBUG(GAout << "EVALUATOR: Nr. of components with min. MSE: " << optNComp + 1 << " (max. " << maxNComp << ")" << std::endl)

		cutoff += trainMSEP.stddev(minNComp) * this->sdfact;

		if(minNComp == 0) {
			optNComp = 1;
		} else {
			optNComp = 0;
			while(optNComp < minNComp && trainMSEP.mean(optNComp) > cutoff) {
				++optNComp;
			}
			if(optNComp <= minNComp) {
				++optNComp;
			}
		}

		IF_DEBUG(GAout << "EVALUATOR: Opt. num. of components: " << optNComp
			<< " (max. " << maxNComp << ")" << std::endl)

		/*
		 * Predict last segment with a model fit to the other observations using optNComp components
		 * (segmentation iterator points to the `test fit` segment
		 */
		this->pls->viewSelectAllRows();
		this->pls->fit(optNComp);

		RSS = arma::accu(arma::square(this->pls->getY() - this->pls->predict(this->pls->getXColumnView(), optNComp)));
	} catch(const ::std::underflow_error& ue) {
		IF_DEBUG(GAout << GAout.lock() << ue.what() << "\n" << GAout.unlock())
		throw Evaluator::EvaluatorException("Can not evaluate variable subset due to an underflow.");
	}

	return RSS;
}

Evaluator* BICEvaluator::clone() const {
	return new BICEvaluator(*this);
}



