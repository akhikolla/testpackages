//
//  Control.h
//
//

#ifndef GenAlgPLS_Control_h
#define GenAlgPLS_Control_h

#include "config.h"

#include <iostream>
#include <vector>
#include <RcppArmadillo.h>

enum VerbosityLevel {
	OFF = 0,
	ON,
	VERBOSE,
	DEBUG_GA,
	DEBUG_EVAL,
	DEBUG_ALL
};

enum CrossoverType {
	SINGLE = 0,
	RANDOM = 1
};

enum FitnessScaling {
	NONE = 0,
	EXP = 1
};

class Control {
public:
	Control(const uint16_t chromosomeSize,
			const uint16_t popSize,
			const uint16_t numGenerations,
			const uint16_t elitism,
			const uint16_t minVariables,
			const uint16_t maxVariables,
			const double mutationProbability,
			const uint16_t numThreads,
			const uint16_t maxDuplicateEliminationTries,
			const double badSolutionThreshold,
			const enum CrossoverType crossover,
			const enum FitnessScaling fitnessScaling,
			const enum VerbosityLevel verbosity) :
	chromosomeSize(chromosomeSize),
	populationSize(popSize),
	numGenerations(numGenerations),
	elitism(elitism),
	minVariables(minVariables),
	maxVariables(maxVariables),
	mutationProbability(mutationProbability),
	numThreads(numThreads),
	maxDuplicateEliminationTries(maxDuplicateEliminationTries),
	badSolutionThreshold(badSolutionThreshold),
	crossover(crossover),
	fitnessScaling(fitnessScaling),
	verbosity(verbosity) {};

	const uint16_t chromosomeSize;
	const uint16_t populationSize;
	const uint16_t numGenerations;
	const uint16_t elitism;
	const uint16_t minVariables;
	const uint16_t maxVariables;
	const double mutationProbability;
	const uint16_t numThreads;
	const uint16_t maxDuplicateEliminationTries;
	const double badSolutionThreshold;
	const enum CrossoverType crossover;
	const enum FitnessScaling fitnessScaling;
	const enum VerbosityLevel verbosity;

	friend std::ostream& operator<<(std::ostream &os, const Control &ctrl) {
		os << "Chromosome size: " << ctrl.chromosomeSize << std::endl
		<< "Population size: " << ctrl.populationSize << std::endl
		<< "Number of generations: " << ctrl.numGenerations << std::endl
		<< "Number of elite chromosomes to keep: " << ctrl.elitism << std::endl
		<< "Number of variables set: " << ctrl.minVariables << " to " << ctrl.maxVariables << std::endl
		<< "Mutation probability: " << ctrl.mutationProbability << std::endl
		<< "Maximum number of tries to eliminate duplicates: " << ctrl.maxDuplicateEliminationTries << std::endl
		<< "Bad solution threshold: " << ctrl.badSolutionThreshold << std::endl
		<< "Crossover-type: " << ((ctrl.crossover == SINGLE) ? "Single" : "Random") << std::endl
		<< "Fitness-scaling: " << ((ctrl.fitnessScaling == EXP) ? "exp" : "None") << std::endl
		<< "Number of threads: " << ctrl.numThreads << std::endl
		<< "Verbosity Level: " << ctrl.verbosity << std::endl
#ifdef ENABLE_DEBUG_VERBOSITY
		<< "Debug enabled"
#else
		<< "Debug disabled"
#endif
		<< std::endl;
		return os;
	};
};

#endif
