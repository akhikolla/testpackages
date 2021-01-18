//
//  MultiThreadedPopulation.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 16.07.2013.
//
//
#include "config.h"

#ifdef HAVE_PTHREAD_H

#include <exception>
#include <vector>
#include <algorithm>
#include <RcppArmadillo.h>
#include <pthread.h>
#include <errno.h>

#include "Logger.h"
#include "RNG.h"
#include "ShuffledSet.h"
#include "MultiThreadedPopulation.h"

using namespace Rcpp;

#ifdef ENABLE_DEBUG_VERBOSITY
#define IF_DEBUG(expr) if(this->ctrl.verbosity == DEBUG_GA || this->ctrl.verbosity == DEBUG_ALL) { expr; }
#define CHECK_PTHREAD_RETURN_CODE(expr) {int rc = expr; if((rc) != 0) { GAerr << "Warning: Call to pthread function failed with error code " << (rc) << " in " << __FILE__ << ":" << __LINE__ << std::endl; }}

#else
#define IF_DEBUG(expr)
#define CHECK_PTHREAD_RETURN_CODE(expr) {expr;}
#endif


/*
 * R user interrupt handling helpers
 */
static inline void check_interrupt_impl(void* /*dummy*/) {
	R_CheckUserInterrupt();
}

inline bool check_interrupt() {
	return (R_ToplevelExec(check_interrupt_impl, NULL) == FALSE);
}

MultiThreadedPopulation::MultiThreadedPopulation(const Control &ctrl, ::Evaluator &evaluator, const std::vector<uint32_t> &seed) : Population(ctrl, evaluator, seed) {
	// initialize original population (generation 0) totally randomly
	if(this->ctrl.numThreads <= 1) {
		throw new std::logic_error("This population should only be used if multiple threads are requested");
	}
	
	this->nextGeneration.reserve(this->ctrl.populationSize);

	int pthreadRC = pthread_mutex_init(&this->syncMutex, NULL);
	if(pthreadRC != 0) {
		throw ThreadingError("Mutex for synchronization could not be initialized");
	}
	
	pthreadRC = pthread_cond_init(&this->startMatingCond, NULL);
	if(pthreadRC != 0) {
		throw ThreadingError("Condition for synchronization (start mating) could not be initialized");
	}
	
	pthreadRC = pthread_cond_init(&this->allThreadsFinishedMatingCond, NULL);
	if(pthreadRC != 0) {
		throw ThreadingError("Condition for synchronization (finished mating) could not be initialized");
	}
	
	this->startMating = false;
	this->allThreadsFinishedMating = false;
	this->killThreads = false;
	
	this->actuallySpawnedThreads = 0;
	this->numThreadsFinishedMating = 0;
}

/**
 * Destructor (destroy mutexes)
 *
 */
MultiThreadedPopulation::~MultiThreadedPopulation() {
	pthread_mutex_destroy(&this->syncMutex);
	pthread_cond_destroy(&this->startMatingCond);
	pthread_cond_destroy(&this->allThreadsFinishedMatingCond);
}

void MultiThreadedPopulation::generateInitialChromosomes(uint16_t numChromosomes, ::Evaluator& evaluator,
		RNG& rng, ShuffledSet& shuffledSet, uint16_t offset, bool checkUserInterrupt) {

	ChVecIt rangeBeginIt = this->nextGeneration.begin() + offset;
	ChVecIt it = rangeBeginIt;
	ChVecIt rangeEndIt = rangeBeginIt + numChromosomes;

	while(it != rangeEndIt && !this->interrupted) {
		(*it) = new Chromosome(this->ctrl, shuffledSet, rng);

		if(std::find_if(rangeBeginIt, it, CompChromsomePtr(*it)) == it) {
			try {
				evaluator.evaluate(**it);
				++it;
			} catch(const ::Evaluator::EvaluatorException &ee) {
				delete (*it);
				if(this->ctrl.verbosity >= VERBOSE) {
					GAout << GAout.lock() << "Could not evaluate chromosome: " << ee.what() << GAout.unlock() << "\n";
				}
			}
		} else {
			delete (*it);
		}

		/*
		 * The main thread has to check for a user interrupt
		 */
		if(checkUserInterrupt == true) {
			GAout.flushThreadSafeBuffer();
			GAerr.flushThreadSafeBuffer();
			if(check_interrupt()) {
				this->interrupted = true;
			}
		}
	}
}

/**
 * Do the actual mating
 *
 */
void MultiThreadedPopulation::mate(uint16_t numChildren, ::Evaluator& evaluator,
	RNG& rng, ShuffledSet& shuffledSet, uint16_t offset,
	bool checkUserInterrupt) {

	double minParentFitness = 0.0;
	
	ChVecIt rangeBeginIt = this->nextGeneration.begin() + offset;
	std::reverse_iterator<ChVecIt> rangeEndIt(rangeBeginIt + numChildren);
	
	Chromosome* tmpChromosome1;
	Chromosome* tmpChromosome2;
	ChVecIt child1It = rangeBeginIt;
	std::reverse_iterator<ChVecIt> child2It = rangeEndIt;
	
	uint8_t child1Tries = 0;
	uint8_t child2Tries = 0;
	std::pair<bool, bool> duplicated(false, false);
	double cutoff = 0.0;
	bool childrenDifferent = true;

	uint32_t discSol1 = 0;
	uint32_t discSol2 = 0;

	while(child1It < child2It.base() && !this->interrupted) {
		childrenDifferent = (child1It + 1 != child2It.base());

		tmpChromosome1 = this->drawChromosomeFromCurrentGeneration(rng(0.0, this->sumCurrentGenFitness));
		do {
			tmpChromosome2 = this->drawChromosomeFromCurrentGeneration(rng(0.0, this->sumCurrentGenFitness));
		} while (tmpChromosome1 == tmpChromosome2);

		tmpChromosome1->mateWith(*tmpChromosome2, rng, *(*child1It), *(*child2It));
		
		minParentFitness = ((tmpChromosome1->getFitness() > tmpChromosome2->getFitness()) ? tmpChromosome1->getFitness() : tmpChromosome2->getFitness());

		(*child1It)->mutate(rng);

		if (childrenDifferent) {
			(*child2It)->mutate(rng);
		}

		/*
		 * Simple rejection
		 * Reject either of the child chromosomes if they are worse than the worst parent (times a given percentage)
		 * or if they are duplicated
		 */
		cutoff = minParentFitness - this->ctrl.badSolutionThreshold * fabs(minParentFitness);

		if(this->ctrl.maxDuplicateEliminationTries > 0) {
			duplicated = Population::checkDuplicated(rangeBeginIt, rangeEndIt, child1It, child2It);
		}

		if((duplicated.first == false) || (++child1Tries > this->ctrl.maxDuplicateEliminationTries)) {
			/*
			 * If the child is a duplicate and we have tried too often
			 * just reset the chromosome to a random point
			 */
			if(duplicated.first == true && child1Tries > this->ctrl.maxDuplicateEliminationTries) {
				(*child1It)->randomlyReset(rng, shuffledSet);
			}

			child1Tries = 0;

			try {
				if(evaluator.evaluate(**child1It) > cutoff) {
					/*
					 * The child is no duplicate (or accepted as one) and is not too bad,
					 * so go on to the next one
					 */
					++child1It;
				} else if(++discSol1 > Population::MAX_DISCARDED_SOLUTIONS_RATIO * numChildren) {
					GAout << GAout.lock() << "Warning: The algorithm may be stuck. Try increasing the badSolutionThreshold!\n" << GAout.unlock();
					discSol1 = 0;
					++child1It;
				}
			} catch(const ::Evaluator::EvaluatorException& ee) {
				if(this->ctrl.verbosity >= VERBOSE) {
					GAout << GAout.lock() << "Could not evaluate chromosome: " << ee.what() << "\n" << GAout.unlock();
				}
			}
		}

		if(childrenDifferent && ((duplicated.second == false) || (++child2Tries > this->ctrl.maxDuplicateEliminationTries))) {
			/*
			 * If the child is a duplicate and we have tried too often
			 * just reset the chromosome to a random point
			 */
			if(duplicated.second == true && child2Tries > this->ctrl.maxDuplicateEliminationTries) {
				(*child2It)->randomlyReset(rng, shuffledSet);
			}

			child2Tries = 0;

			try {
				if(evaluator.evaluate(**child2It) > cutoff) {
					/*
					 * The child is no duplicate (or accepted as one) and is not too bad,
					 * so go on to the next one
					 */
					++child2It;
				} else if(++discSol2 > Population::MAX_DISCARDED_SOLUTIONS_RATIO * numChildren) {
					GAout << GAout.lock() << "Warning: The algorithm may be stuck. Try increasing the badSolutionThreshold!\n" << GAout.unlock();
					discSol2 = 0;
					++child2It;
				}
			} catch(const ::Evaluator::EvaluatorException& ee) {
				if(this->ctrl.verbosity >= VERBOSE) {
					GAout << GAout.lock() << "Could not evaluate chromosome: " << ee.what() << "\n" << GAout.unlock();
				}
			}
		}

		/*
		 * The main thread has to check for a user interrupt
		 */
		if(checkUserInterrupt == true) {
			GAout.flushThreadSafeBuffer();
			GAerr.flushThreadSafeBuffer();
			if(check_interrupt()) {
				this->interrupted = true;
			}
		}
	}
}

/**
 * Start the evolution
 */
void MultiThreadedPopulation::run() {
	int i = 0, j = 0;
	RNG rng(this->seed);
	double minFitness = 0.0;
	ShuffledSet shuffledSet(this->ctrl.chromosomeSize);
	MultiThreadedPopulation::ThreadArgsWrapper* threadArgs;
	uint16_t maxThreadsToSpawn = this->ctrl.numThreads - 1;
	uint16_t numChildrenPerThread = this->ctrl.populationSize / this->ctrl.numThreads;
	int remainingChildren = this->ctrl.populationSize % this->ctrl.numThreads;
	uint16_t numChildrenMainThread = numChildrenPerThread;
	uint16_t offset = 0;
	pthread_attr_t threadAttr;
	pthread_t* threads;

	/*****************************************************************************************
	 * Initialize the current/next generation and enable thread safety for the output
	 *****************************************************************************************/
	if(this->ctrl.verbosity > OFF) {
		GAout << "Generating initial population" << std::endl;
	}

	/* 1. fill the next generation with null pointers and initialize the current generation */
	this->nextGeneration.resize(this->ctrl.populationSize, (Chromosome*) NULL);
	this->initCurrentGeneration(shuffledSet, rng);


	/* 2. let the threads generate and evaluate a bunch of chromosomes ... */
	GAout.enableThreadSafety(true);
	GAerr.enableThreadSafety(true);

	/*****************************************************************************************
	 * Setup threads
	 *****************************************************************************************/
	
	threadArgs = new MultiThreadedPopulation::ThreadArgsWrapper[maxThreadsToSpawn];
	threads = new pthread_t[maxThreadsToSpawn];
	
	int pthreadRC = pthread_attr_init(&threadAttr);
	if(pthreadRC != 0) {
		throw ThreadingError("Thread attributes could not be initialized");
	}
	
	pthreadRC = pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);
	
	if(pthreadRC != 0) {
		throw ThreadingError("Thread attributes could not be modified to make the thread joinable");
	}
	
	for(i = maxThreadsToSpawn - 1; i >= 0; --i) {
		threadArgs[i].numChildren = numChildrenPerThread;
		
		if(remainingChildren > 0) {
			--remainingChildren;
			++threadArgs[i].numChildren;
		}
		threadArgs[i].offset = offset;
		threadArgs[i].popObj = this;
		threadArgs[i].seed = rng();
		threadArgs[i].evalObj = this->evaluator.clone();
		threadArgs[i].chromosomeSize = this->ctrl.chromosomeSize;

		/*
		 * Once created, the threads already start generating the initial generation!
		 */
		pthreadRC = pthread_create((threads + i), &threadAttr, &MultiThreadedPopulation::matingThreadStart, (void *) (threadArgs + i));
		
		if(pthreadRC == 0) {
			++this->actuallySpawnedThreads;
			offset += threadArgs[i].numChildren;
		} else {
			numChildrenMainThread += threadArgs[i].numChildren;
			IF_DEBUG(GAerr << "Warning: Thread " << i << " could not be created: " << strerror(pthreadRC) << std::endl;)
		}
	}
	
	CHECK_PTHREAD_RETURN_CODE(pthread_attr_destroy(&threadAttr))
	
	if(this->actuallySpawnedThreads < maxThreadsToSpawn) {
		GAerr << GAerr.lock() << "Warning: Only " << this->actuallySpawnedThreads << " threads could be spawned\n" << GAerr.unlock();
	} else if(this->ctrl.verbosity >= ON) {
		GAout  << GAout.lock() << "Spawned " << this->actuallySpawnedThreads << " threads\n" << GAout.unlock();
	}

	/*****************************************************************************************
	 * Generate initial population
	 *****************************************************************************************/

	this->generateInitialChromosomes(numChildrenMainThread, this->evaluator, rng, shuffledSet, offset, true);

	/*****************************************************************************************
	 * Wait for threads to create current generation and further process the initial population
	 *****************************************************************************************/
	this->waitForAllThreadsToFinishMating();

	/* Maybe check the initial generation for duplicats ??? */

	/*
	 * Signal output streams that multithreading is over
	 */
	GAout.enableThreadSafety(false);
	GAerr.enableThreadSafety(false);

	if(this->interrupted == false) {
		/***********************************************************************
		 * Update minFitness, the current generation and the sumFitness
		 * and print the generation if requested
		 **********************************************************************/
		minFitness = (*(std::min_element(this->nextGeneration.begin(), this->nextGeneration.end(), MultiThreadedPopulation::OrderChromosomePtr())))->getFitness();

		this->sumCurrentGenFitness = this->updateCurrentGeneration(this->nextGeneration, minFitness, true, true);

		if(this->ctrl.verbosity >= VERBOSE && this->ctrl.verbosity != DEBUG_EVAL) {
			this->printCurrentGeneration();
		}
	}

	/*****************************************************************************************
	 * Generate remaining generations
	 *****************************************************************************************/
	
	for(i = this->ctrl.numGenerations; i > 0 && !this->interrupted; --i) {
		IF_DEBUG(GAout << "Unique chromosomes: " << this->countUniques() << std::endl;)
		
		if(this->ctrl.verbosity > OFF) {
			GAout << "Generating generation " << (this->ctrl.numGenerations - i + 1) << std::endl;
		}

		/*
		 * Prepare output streams for multithreading again
		 */
		GAout.enableThreadSafety(true);
		GAerr.enableThreadSafety(true);

		/*****************************************************************************************
		 * broadcast to all threads to start mating
		 *****************************************************************************************/
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->syncMutex))
		
		this->startMating = true;
		
		CHECK_PTHREAD_RETURN_CODE(pthread_cond_broadcast(&this->startMatingCond))
		
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->syncMutex))
		
		/*
		 * Mate two chromosomes to generate two children that are eventually mutated
		 * To get the same population size, a total of popSize / 2 mating pairs have
		 * to generate 2 children
		 *
		 */
		this->mate(numChildrenMainThread, this->evaluator, rng, shuffledSet, offset, true);
		
		this->waitForAllThreadsToFinishMating();
		/*
		 * Signal output streams that multithreading is over
		 */
		GAout.enableThreadSafety(false);
		GAerr.enableThreadSafety(false);

		/***********************************************************************
		 * Update minFitness, the current generation and the sumFitness
		 * and print the generation if requested
		 **********************************************************************/
		minFitness = (*(std::min_element(this->nextGeneration.begin(), this->nextGeneration.end(), MultiThreadedPopulation::OrderChromosomePtr())))->getFitness();

		this->sumCurrentGenFitness = this->updateCurrentGeneration(this->nextGeneration, minFitness, false, true);

		if(this->ctrl.verbosity >= VERBOSE && this->ctrl.verbosity != DEBUG_EVAL) {
			this->printCurrentGeneration();
		}
	}
	
	/*****************************************************************************************
	 * Signal threads to end
	 *****************************************************************************************/
	GAout.enableThreadSafety(true);
	GAerr.enableThreadSafety(true);
	CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->syncMutex))
	
	this->startMating = true;
	this->killThreads = true;
	
	CHECK_PTHREAD_RETURN_CODE(pthread_cond_broadcast(&this->startMatingCond))
	
	CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->syncMutex))
	
	
	for(i = maxThreadsToSpawn - 1; i >= 0; --i) {
		/*
		 * If the thread was never created (i.e. pthread_create failed) the call will return an
		 * error code, but it will not block the thread!
		 */
		CHECK_PTHREAD_RETURN_CODE(pthread_join(threads[i], NULL))
		
		delete threadArgs[i].evalObj;
	}
	
	delete threadArgs;
	delete threads;

	GAout.enableThreadSafety(false);
	GAerr.enableThreadSafety(false);

	/*****************************************************************************************
	 * Update elite and the generation
	 * and delete the old `nextGeneration`
	 *****************************************************************************************/
	
	for(j = 0; j < this->ctrl.populationSize; ++j) {
		if(this->nextGeneration[j]) {
			delete this->nextGeneration[j];
		}
	}
}


/**
 * Setup and start the mating threads
 */
void* MultiThreadedPopulation::matingThreadStart(void* obj) {
	ThreadArgsWrapper* args = static_cast<ThreadArgsWrapper*>(obj);
	RNG rng(args->seed);
	ShuffledSet shuffledSet(args->chromosomeSize);

	/* First generate a bunch of initial chromosomes */
	args->popObj->generateInitialChromosomes(args->numChildren, *args->evalObj, rng, shuffledSet, args->offset, false);
	args->popObj->waitForAllThreadsToFinishMating();

	/* The start the mating cycle */
	args->popObj->runMating(args->numChildren, *args->evalObj, rng, shuffledSet, args->offset);
	return NULL;
}

/**
 * Run the mating control loop
 */
void MultiThreadedPopulation::runMating(uint16_t numMatingCouples, ::Evaluator& evaluator,
		RNG& rng, ShuffledSet& shuffledSet, uint16_t offset) {
	while(true) {
		/*****************************************************************************************
		 * Wait until the thread is started
		 *****************************************************************************************/
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->syncMutex))
		
		while(this->startMating == false) {
			CHECK_PTHREAD_RETURN_CODE(pthread_cond_wait(&this->startMatingCond, &this->syncMutex))
		}
		
		/*****************************************************************************************
		 * Check if the thread is killed
		 *****************************************************************************************/
		if(this->killThreads == true) {
			CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->syncMutex))
			break;
		}
		
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->syncMutex))
		
		/*****************************************************************************************
		 * Do actual mating
		 *****************************************************************************************/
		this->mate(numMatingCouples, evaluator, rng, shuffledSet, offset, false);
		
		/*****************************************************************************************
		 * Signal that the thread has finished mating
		 *****************************************************************************************/
		this->waitForAllThreadsToFinishMating();
	}
}

/**
 * Wait for all threads to finish the current generation
 */
inline void MultiThreadedPopulation::waitForAllThreadsToFinishMating() {
	CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->syncMutex))
	
	if(++this->numThreadsFinishedMating > this->actuallySpawnedThreads) { // > because the main thread must finish mating as well
		this->allThreadsFinishedMating = true;
		this->numThreadsFinishedMating = 0;
		this->startMating = false;
		
		CHECK_PTHREAD_RETURN_CODE(pthread_cond_broadcast(&this->allThreadsFinishedMatingCond))
	} else {
		this->allThreadsFinishedMating = false;
	}
	
	//	pthreadRC = pthread_mutex_unlock(&this->syncMutex);
	//	CHECK_PTHREAD_RETURN_CODE(pthreadRC)
	//	int pthreadRC = pthread_mutex_lock(&this->syncMutex);
	//	CHECK_PTHREAD_RETURN_CODE(pthreadRC)
	
	while(this->allThreadsFinishedMating == false) {
		CHECK_PTHREAD_RETURN_CODE(pthread_cond_wait(&this->allThreadsFinishedMatingCond, &this->syncMutex))
	}
	
	CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->syncMutex))
}

#endif
