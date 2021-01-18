//
//  GenAlg.cpp
//  GenAlgPLS
//

#include "config.h"

#include <RcppArmadillo.h>
#include <set>
#include <memory>

#include "Logger.h"
#include "Chromosome.h"
#include "Control.h"
#include "PLS.h"
#include "UserFunEvaluator.h"
#include "PLSEvaluator.h"
#include "LMEvaluator.h"
#include "BICEvaluator.h"
#include "SingleThreadPopulation.h"
#include "RNG.h"

#ifdef HAVE_PTHREAD_H
#include "MultiThreadedPopulation.h"
#endif

#include "GenAlg.h"

using namespace Rcpp;

/**
 * .C entry point definitions for R
 */
static const R_CallMethodDef exportedCallMethods[] = {
    {"C_genAlgPLS", (DL_FUNC) &genAlgPLS, 4},
    {"C_evaluate", (DL_FUNC) &evaluate, 5},
    {"C_simpls", (DL_FUNC) &simpls, 5},
    {NULL, NULL, 0}
};


extern "C" void R_init_gaselect(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, exportedCallMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

RcppExport SEXP genAlgPLS(SEXP Scontrol, SEXP SX, SEXP Sy, SEXP Sseed) {
  std::unique_ptr<::Evaluator> eval;
	std::unique_ptr<PLS> pls;
	std::unique_ptr<Population> pop;
BEGIN_RCPP
	List control = List(Scontrol);
	uint32_t singleSeed = as<uint32_t>(Sseed);
	std::vector<uint32_t> seed;
	uint16_t numThreads = as<uint16_t>(control["numThreads"]);
	VerbosityLevel verbosity = (VerbosityLevel) as<int>(control["verbosity"]);
	EvaluatorClass evalClass = (EvaluatorClass) as<int>(control["evaluatorClass"]);

#ifdef ENABLE_DEBUG_VERBOSITY
	PLSEvaluator::counter = 0;
#endif

	if(numThreads > 1) {
#ifdef HAVE_PTHREAD_H
		if(evalClass == USER) {
			GAerr << "Warning: Multithreading is not available when using a user supplied function for evaluation" << std::endl;
		}
#else
		GAerr << "Warning: Threads are not supported on this system" << std::endl;
		numThreads = 1;
#endif
	} else if(numThreads < 1) {
		numThreads = 1;
	}

	// All checks are disabled and must be performed in the R code calling this script
	// Otherwise unexpected behaviour
	Control ctrl(as<uint16_t>(control["chromosomeSize"]),
				 as<uint16_t>(control["populationSize"]),
				 as<uint16_t>(control["numGenerations"]),
				 as<uint16_t>(control["elitism"]),
				 as<uint16_t>(control["minVariables"]),
				 as<uint16_t>(control["maxVariables"]),
				 as<double>(control["mutationProb"]),
				 numThreads,
				 as<uint16_t>(control["maxDuplicateEliminationTries"]),
				 as<double>(control["badSolutionThreshold"]),
				 (CrossoverType) as<int>(control["crossover"]),
				 (FitnessScaling) as<int>(control["fitnessScaling"]),
				 verbosity);

	/*
	 * Generate a common seed for the Population and the PLSEvaluator objects
	 */
	RNG rng(singleSeed);
	seed.reserve(RNG::SEED_SIZE);
	for(uint32_t i = 0; i < RNG::SEED_SIZE; ++i) {
		seed.push_back(rng());
	}

	switch(evalClass) {
		case USER: {
			eval.reset(new UserFunEvaluator(as<Rcpp::Function>(control["userEvalFunction"]), ctrl.verbosity));
			break;
		}
		case PLS_EVAL: {
			Rcpp::NumericMatrix XMat(SX);
			Rcpp::NumericMatrix YMat(Sy);
			arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
			arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
			PLSMethod method = (PLSMethod) as<int>(control["plsMethod"]);

			pls = PLS::getInstance(method, X, Y.col(0));

			eval.reset(new PLSEvaluator(std::move(pls), as<uint16_t>(control["numReplications"]),
                               as<uint16_t>(control["maxNComp"]), seed, ctrl.verbosity,
                               as<uint16_t>(control["innerSegments"]),
                               as<uint16_t>(control["outerSegments"]),
                               as<double>(control["testSetSize"]),
                               as<double>(control["sdfact"])));

			break;
		}
		case PLS_FIT: {
			Rcpp::NumericMatrix XMat(SX);
			Rcpp::NumericMatrix YMat(Sy);
			arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
			arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
			PLSMethod method = (PLSMethod) as<int>(control["plsMethod"]);

			pls = PLS::getInstance(method, X, Y.col(0));

			BICEvaluator::Statistic stat = (BICEvaluator::Statistic) as<int>(control["statistic"]);

			eval.reset(new BICEvaluator(std::move(pls),
                               as<uint16_t>(control["maxNComp"]),
                               seed,
                               ctrl.verbosity,
                               as<uint16_t>(control["innerSegments"]),
                               stat,
                               as<double>(control["sdfact"])));

			break;
		}
		case LM: {
			Rcpp::NumericMatrix XMat(SX);
			Rcpp::NumericMatrix YMat(Sy);
			arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
			arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
			arma::colvec y = Y.col(0);

			LMEvaluator::Statistic stat = (LMEvaluator::Statistic) as<int>(control["statistic"]);
			eval.reset(new LMEvaluator(X, y, stat, verbosity));

			break;
		}
		default:
			throw Rcpp::exception("No valid evaluation method was selected.", __FILE__, __LINE__);
			break;
	}

	if(ctrl.verbosity >= VERBOSE) {
		GAout << ctrl << std::endl;
	}

#ifdef HAVE_PTHREAD_H
	try {
		if(numThreads > 1) {
			pop.reset(new MultiThreadedPopulation(ctrl, *eval, seed));
		} else {
			pop.reset(new SingleThreadPopulation(ctrl, *eval, seed));
		}
		pop->run();
	} catch(MultiThreadedPopulation::ThreadingError& te) {
		if(ctrl.verbosity >= DEBUG_GA) {
			throw te;
		} else {
			throw Rcpp::exception("Multithreading could not be initialized. Set numThreads to 0 to avoid this problem.", __FILE__, __LINE__);
		}
	}
#else
	pop.reset(new SingleThreadPopulation(ctrl, *eval, seed));
	pop->run();
#endif

#ifdef ENABLE_DEBUG_VERBOSITY
	Rcpp::Rcout << "Called evaluator " << PLSEvaluator::counter << " times" << std::endl;
#endif

	if(pop->wasInterrupted() == true) {
		GAout << "Interrupted - returning best solutions found so far" << std::endl;
	}

	/*
	 * Gather data for return
	 */

	Population::SortedChromosomes result = pop->getResult();

	Rcpp::NumericVector retFitnessEvolution(pop->getFitnessEvolution().begin(), pop->getFitnessEvolution().end());
	std::vector<arma::uvec> segmentation = eval->getSegmentation();
	Rcpp::LogicalMatrix retMatrix(ctrl.chromosomeSize, (const int) result.size());
	Rcpp::NumericVector retFitnesses((const int) result.size());
	uint16_t i = (uint16_t) result.size() - 1;

	for(Population::SortedChromosomes::iterator it = result.begin(); it != result.end(); ++it, --i) {
		retFitnesses[i] = it->getFitness();
		retMatrix.column(i) = it->toLogicalVector();
	}


	return Rcpp::List::create(Rcpp::Named("subsets") = retMatrix,
							  Rcpp::Named("fitness") = retFitnesses,
							  Rcpp::Named("fitnessEvolution") = retFitnessEvolution,
							  Rcpp::Named("segmentation") = Rcpp::wrap(segmentation));
VOID_END_RCPP
	return R_NilValue;
}

/**
 *
 */
RcppExport SEXP evaluate(SEXP Sevaluator, SEXP SX, SEXP Sy, SEXP Ssubsets, SEXP Sseed) {
	std::unique_ptr<::Evaluator> eval;
  std::unique_ptr<PLS> pls;

	BEGIN_RCPP
	List evaluator = List(Sevaluator);
	Rcpp::NumericMatrix XMat(SX);
	Rcpp::NumericMatrix YMat(Sy);
	Rcpp::LogicalMatrix subsets(Ssubsets);
	Rcpp::NumericVector fitness(subsets.cols());
	std::vector<arma::uvec> segmentation;
	arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
	arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
	EvaluatorClass evalClass = (EvaluatorClass) as<int>(evaluator["evaluatorClass"]);
	std::vector<uint32_t> seed;

	switch(evalClass) {
		case USER: {
			eval.reset(new UserFunEvaluator(as<Rcpp::Function>(evaluator["userEvalFunction"]), OFF));
			break;
		}
		case PLS_EVAL: {
			PLSMethod method = (PLSMethod) as<int>(evaluator["plsMethod"]);
			RNG rng(as<uint32_t>(Sseed));

			seed.reserve(RNG::SEED_SIZE);
			for (uint32_t i = 0; i < RNG::SEED_SIZE; ++i) {
				seed.push_back(rng());
			}

			pls = PLS::getInstance(method, X, Y.col(0));

			eval.reset(new PLSEvaluator(std::move(pls), as<uint16_t>(evaluator["numReplications"]),
                               as<uint16_t>(evaluator["maxNComp"]), seed,
                               (VerbosityLevel) as<int>(evaluator["verbosity"]),
                               as<uint16_t>(evaluator["innerSegments"]),
                               as<uint16_t>(evaluator["outerSegments"]),
                               as<double>(evaluator["testSetSize"]),
                               as<double>(evaluator["sdfact"])));

			break;
		}
		case PLS_FIT: {
			Rcpp::NumericMatrix XMat(SX);
			Rcpp::NumericMatrix YMat(Sy);
			arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
			arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
			PLSMethod method = (PLSMethod) as<int>(evaluator["plsMethod"]);
			RNG rng(as<uint32_t>(Sseed));

			seed.reserve(RNG::SEED_SIZE);
			for (uint32_t i = 0; i < RNG::SEED_SIZE; ++i) {
				seed.push_back(rng());
			}

			pls = PLS::getInstance(method, X, Y.col(0));

			BICEvaluator::Statistic stat = (BICEvaluator::Statistic) as<int>(evaluator["statistic"]);

			eval.reset(new BICEvaluator(std::move(pls),
                               as<uint16_t>(evaluator["maxNComp"]),
                               seed,
                               (VerbosityLevel) as<int>(evaluator["verbosity"]),
                               as<uint16_t>(evaluator["innerSegments"]),
                               stat,
                               as<double>(evaluator["sdfact"])));

			break;
		}
		case LM: {
			arma::colvec y = Y.col(0);

			LMEvaluator::Statistic stat = (LMEvaluator::Statistic) as<int>(evaluator["statistic"]);
			eval.reset(new LMEvaluator(X, y, stat, (VerbosityLevel) as<int>(evaluator["verbosity"])));
			break;
		}
		default:
			break;
	}
	int row = 0;
	uint16_t i = 0;

	for(int col = 0; col < subsets.cols(); ++col) {
		arma::uvec selectedColumns(subsets.rows());
		i = 0;

		for(row = 0; row < subsets.rows(); ++row) {
			if(subsets(row, col) == true) {
				selectedColumns[i++] = row;
			}
		}
		if(i > 0) {
			selectedColumns.resize(i);
			fitness[col] = eval->evaluate(selectedColumns);
		}
	}

	segmentation = eval->getSegmentation();

	return Rcpp::List::create(Rcpp::Named("fitness") = Rcpp::wrap(fitness),
							  Rcpp::Named("segmentation") = Rcpp::wrap(segmentation));

	VOID_END_RCPP
	return R_NilValue;
}

//RcppExport SEXP WELL19937a(SEXP Sn, SEXP Sseed) {
//	uint32_t n = as<uint32_t>(Sn);
//	uint32_t seed = as<uint32_t>(Sseed);
//
//	RNG rng(seed);
//
//	Rcpp::NumericVector retMat(n);
//
//	for(uint16_t row = 0; row < n; ++row) {
//		retMat[row] = rng();
//	}
//
//	return Rcpp::wrap(retMat);
//}

RcppExport SEXP simpls(SEXP Xs, SEXP Ys, SEXP ncomps, SEXP newXs, SEXP reps) {
BEGIN_RCPP
	Rcpp::NumericMatrix XMat(Xs);
	Rcpp::NumericVector YVec(Ys);
	Rcpp::NumericMatrix newXMat(newXs);
	uint16_t ncomp = Rcpp::as<uint16_t>(ncomps);
	int rep = Rcpp::as<int>(reps);

	arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
	arma::vec Y(YVec.begin(), YVec.length(), false);
	arma::mat newX(newXMat.begin(), newXMat.nrow(), newXMat.ncol(), false);

	PLSSimpls simpls(X, Y);
	for(; rep >= 0; --rep) {
		simpls.fit(ncomp);
	}

	return Rcpp::List::create(Rcpp::Named("coefficients") = simpls.getCoefficients(),
							  Rcpp::Named("predicted") = simpls.predict(newX));
END_RCPP
}

//
//
// RcppExport SEXP evalTest(SEXP Xs, SEXP Ys, SEXP numReplications, SEXP numSegments) {
// BEGIN_RCPP
// 	Rcpp::NumericMatrix XMat(Xs);
// 	Rcpp::NumericMatrix YMat(Ys);
// 	arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
// 	arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
// 	PLSSimpls pls(X, Y, false);
// 	PLSEvaluator eval(pls, Rcpp::as<uint16_t>(numReplications), Rcpp::as<uint16_t>(numSegments), DEBUG_VERBOSE);
//
// 	arma::uvec colSubset(X.n_cols);
//
// 	for(uint16_t i = 0; i < X.n_cols; ++i) {
// 		colSubset[i] = i;
// 	}
//
// 	return Rcpp::wrap(eval.evaluate(colSubset));
// END_RCPP
// }
//
