//
//  GenAlg.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 05.04.2013.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef GenAlgPLS_GenAlg_h
#define GenAlgPLS_GenAlg_h

#include "config.h"

#include <RcppArmadillo.h>

#ifdef TIMING_BENCHMARK

#include <sys/time.h>
#include <sys/types.h>

#endif

enum EvaluatorClass {
	USER = 0,
	PLS_EVAL = 1,
	LM = 2,
	PLS_FIT = 3
};

/**
 * arguments:
 *	control ... A R list with following entries:
 *		uint16_t chromosomeSize ... The size of the chromosome (most times equal to the number of columns of X) (> 0)
 *		uint16_t populationSize ... Number of indivudual chromosomes in the population (i.e. per generation) (> 0)
 *		uint16_t numGenerations ... The number of generations to generate (> 0)
 *		uint16_t minVariables ... The minimum number of variables in a subset
 *		uint16_t maxVariables ... The maximum number of variables in a subset
 *		uint16_t elitism ... The number of "elite" chromosomes to keep accross all generations (>= 0)
 *		double mutationProb ... The probability of using a new variable (0 <= onesRatio < 1)
 *		uint16_t numThreads ... The maximum number of threads to spawn
 *		uint16_t maxDuplicateEliminationTries ... The maximum number of tries to eliminate duplicates
 *		double badSolutionThreshold ... The better child must not be more than badSolutionThreadshold percent worse than the worse parent
 *		CrossoverType crossover ... Type of crossover to use
 *		FitnessScaling fitnessScaling ... How to scale the fitness (0 = NONE, 1 = EXP)
 *		VerbosityLevel verbosity ... Level of verbosity
 *		EvaluatorClass evaluatorClass ... The evaluator to use
 *		Rcpp::Function userEvalFunction ... The function to be called for evaluating the fitness of a chromosome
 *		PLSMethod plsMethod ... PLS method to use in internal evaluation
 *		uint16_t numReplications ... Number of replications in the internal evaluation procedure
 *			(the variable subset is evaluted with CV numReplication times and the mean fitness is returned) (> 0)
 *		uint16_t innerSegments ... Number of CV segments used in the internal evaluation method (> 0)
 *		uint16_t outerSegments ... Number of outer CV segments (0 or 1 = use srCV, > 0 use the rdCV strategy)
 *		double testSetSize ... If srCV should be used, the rel. size of the test set between 0 and 1 (ignored if outerSegments > 1)
 *		double sdfact ... The factor to scale the SD with when selecting the optimal number of components
 *		uint16_t maxNComp ... The maximum number of componentes the PLS models should consider
 *		int statistic ... The statistic the LM Evaluator should use
 *	X ... A numeric matrix with dimensions n x p (optional - only needed if using internal evaluation methods)
 *	y ... A numeric vector with length n (optional - only needed if using internal evaluation methods)
 *	seed ... An integer (uint32_t) with the initial seed
 */
RcppExport SEXP genAlgPLS(SEXP control, SEXP X, SEXP y, SEXP seed);

/**
 * evaluate the given data with the given evaluator
 * arguments:
 *	evaluator ... A R list with following entries
 *		VerbosityLevel verbosity ... Level of verbosity
 *		EvaluatorClass evaluatorClass ... The evaluator to use
 *		Rcpp::Function userEvalFunction ... The function to be called for evaluating the fitness of a chromosome
 *		PLSMethod plsMethod ... PLS method to use in internal evaluation
 *		uint16_t numReplications ... Number of replications in the internal evaluation procedure
 *			(the variable subset is evaluted with CV numReplication times and the mean fitness is returned) (> 0)
 *		uint16_t innerSegments ... Number of CV segments used in the internal evaluation method (> 0)
 *		uint16_t outerSegments ... Number of outer CV segments (0 or 1 = use srCV, > 0 use the rdCV strategy)
 *		double testSetSize ... If srCV should be used, the rel. size of the test set between 0 and 1 (ignored if outerSegments > 1)
 *		double sdfact ... The factor to scale the SD with when selecting the optimal number of components
 *		uint16_t maxNComp ... The maximum number of componentes the PLS models should consider
 *		int sepTransformation ... The type of transformation for the SEP (NONE or LOG)
 *		int statistic ... The statistic the LM Evaluator should use
 *
 *	X ... A numeric matrix with dimensions n x p
 *	y ... A numeric vector with length n
 *	subsets ... A logical matrix with dimensions p x k, where k is the number of different subsets to evaluate
 *	seed ... An integer (uint32_t) with the initial seed
 */
RcppExport SEXP evaluate(SEXP evaluator, SEXP X, SEXP y, SEXP subsets, SEXP seed);

RcppExport SEXP simpls(SEXP X, SEXP y, SEXP ncomp, SEXP newX, SEXP rep);
//RcppExport SEXP WELL19937a(SEXP n, SEXP seed);

#endif
