//
//  ThreadedPopulation.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 16.07.2013.
//
//

#ifndef GenAlgPLS_ThreadedPopulation_h
#define GenAlgPLS_ThreadedPopulation_h
#include "config.h"

#ifdef HAVE_PTHREAD_H

#include <stdexcept>
#include <iostream>
#include <streambuf>
#include <vector>
#include <set>
#include <string>
#include <utility>

#include "Chromosome.h"
#include "Evaluator.h"
#include "Control.h"
#include "Population.h"

#include "RNG.h"

class MultiThreadedPopulation : public Population {
	
public:
	MultiThreadedPopulation(const Control &ctrl, ::Evaluator &evaluator, const std::vector<uint32_t> &seed);
	~MultiThreadedPopulation();
	
	/**
	 * Returns the last generation and the elite (if any)
	 * Invalid once the Population object is destroyed!
	 */
	
	class ThreadingError : public std::runtime_error {
	public:
		ThreadingError(const char* what) : std::runtime_error(what) {};
		virtual ~ThreadingError() throw() {};
	};
	
	void run();
private:
	struct ThreadArgsWrapper {
		MultiThreadedPopulation* popObj;
		Evaluator* evalObj;
		uint32_t seed;
		uint16_t numChildren;
		uint16_t offset;
		uint16_t chromosomeSize;
	};

	ChVec nextGeneration;
	double sumCurrentGenFitness;
	
	/*
	 * Mutex and condition variables
	 */
	pthread_mutex_t syncMutex;
	pthread_cond_t startMatingCond;
	pthread_cond_t allThreadsFinishedMatingCond;
	
	bool startMating;
	bool killThreads;
	bool allThreadsFinishedMating;

	uint16_t actuallySpawnedThreads;
	uint16_t numThreadsFinishedMating;

	inline void generateInitialChromosomes(uint16_t numChromosomes, ::Evaluator& evaluator,
		RNG& rng, ShuffledSet& shuffledSet, uint16_t offset,
		bool checkUserInterrupt = true);

	inline void mate(uint16_t numChildren, ::Evaluator& evaluator,
		RNG& rng, ShuffledSet& shuffledSet, uint16_t offset,
		bool checkUserInterrupt = true);
	
	static void* matingThreadStart(void* obj);

	inline void runMating(uint16_t numMatingCoupls, ::Evaluator& evaluator,
		RNG& rng, ShuffledSet& shuffledSet, uint16_t offset);

	inline void waitForAllThreadsToFinishMating();
	
	class OrderChromosomePtr : public std::binary_function<Chromosome*, Chromosome*, bool> {
	public:
		bool operator()(const Chromosome* const ch1, const Chromosome* const ch2) {
			return (ch1->getFitness() < ch2->getFitness());
		};
	};
};


#endif
#endif
