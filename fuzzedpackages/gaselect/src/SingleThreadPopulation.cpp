//
//  SingleThreadPopulation.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 16.07.2013.
//
//
#include "config.h"

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <RcppArmadillo.h>

#include "Logger.h"
#include "RNG.h"
#include "SingleThreadPopulation.h"
#include "ShuffledSet.h"

using namespace Rcpp;

#ifdef ENABLE_DEBUG_VERBOSITY
#define IF_DEBUG(expr) if(this->ctrl.verbosity == DEBUG_GA || this->ctrl.verbosity == DEBUG_ALL) { expr; }
#else
#define IF_DEBUG(expr)
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

SingleThreadPopulation::SingleThreadPopulation(const Control &ctrl, ::Evaluator &evaluator, const std::vector<uint32_t> &seed) : Population(ctrl, evaluator, seed) {}


void SingleThreadPopulation::run() {
	int i = 0;
	ShuffledSet shuffledSet(this->ctrl.chromosomeSize);
	RNG rng(this->seed);
	
	double sumFitness = 0.0;
	double minFitness = 0.0;
	double minParentFitness = 0.0;
	
	ChVec newGeneration;
	
	Chromosome* tmpChromosome1;
	Chromosome* tmpChromosome2;
	ChVecIt child1It;
	ChVecRIt child2It;
	uint8_t child1Tries = 0;
	uint8_t child2Tries = 0;
	std::pair<bool, bool> duplicated(false, false);
	double cutoff = 0.0;
	bool childrenDifferent = true;

	uint32_t discSol1 = 0;
	uint32_t discSol2 = 0;
	uint32_t maxDiscardedSolutions = Population::MAX_DISCARDED_SOLUTIONS_RATIO * this->ctrl.populationSize;

	newGeneration.reserve(this->ctrl.populationSize);
	
	if(this->ctrl.verbosity > OFF) {
		GAout << "Generating initial population" << std::endl;
	}
	
	while(newGeneration.size() < this->ctrl.populationSize && !this->interrupted) {
		tmpChromosome1 = new Chromosome(this->ctrl, shuffledSet, rng);
		
		/* Check if chromosome is already in the initial population */
		if(std::find_if(newGeneration.begin(), newGeneration.end(), CompChromsomePtr(tmpChromosome1)) == newGeneration.end()) {
			try {
				this->evaluator.evaluate(*tmpChromosome1);

				if(tmpChromosome1->getFitness() < minFitness) {
					minFitness = tmpChromosome1->getFitness();
				}
				
				this->addChromosomeToElite(*tmpChromosome1);
				
				newGeneration.push_back(tmpChromosome1);
			} catch(const ::Evaluator::EvaluatorException& ee) {
				delete tmpChromosome1;
				if(this->ctrl.verbosity >= VERBOSE) {
					GAout << "Could not evaluate chromosome: " << ee.what() << std::endl;
				}
			}
		} else {
			delete tmpChromosome1;
		}

		if(check_interrupt()) {
			this->interrupted = true;
		}
	}

	this->initCurrentGeneration(shuffledSet, rng);

	/*
	 * Transform the fitness map of the current generation to start at 0
	 * and copy old generation to new generation
	 */
	sumFitness = this->updateCurrentGeneration(newGeneration, minFitness, true);

	if(this->ctrl.verbosity >= VERBOSE && this->ctrl.verbosity != DEBUG_EVAL) {
		this->printCurrentGeneration();
	}
	
	for(i = this->ctrl.numGenerations; i > 0 && !this->interrupted; --i) {
		minFitness = 0.0;

		IF_DEBUG(
			GAout << "Unique chromosomes: " << this->countUniques() << std::endl;
		)

		if(this->ctrl.verbosity > OFF) {
			GAout << "Generating generation " << (this->ctrl.numGenerations - i + 1) << std::endl;
		}
		
		child1It = newGeneration.begin();
		child2It = newGeneration.rbegin();



		while(child1It < child2It.base() && !this->interrupted) {
			childrenDifferent = (child1It + 1 != child2It.base());

			tmpChromosome1 = this->drawChromosomeFromCurrentGeneration(rng(0.0, sumFitness));
			do {
				tmpChromosome2 = this->drawChromosomeFromCurrentGeneration(rng(0.0, sumFitness));
			} while(tmpChromosome1 == tmpChromosome2);

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
				duplicated = Population::checkDuplicated(newGeneration.begin(), newGeneration.rbegin(), child1It, child2It);
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
					if((this->evaluator.evaluate(**child1It) > cutoff) || (discSol1 > maxDiscardedSolutions)) {
						if((*child1It)->getFitness() < minFitness) {
							minFitness = (*child1It)->getFitness();
						}

						this->addChromosomeToElite(**child1It);

						/*
						 * The child is no duplicate (or accepted as one) and is not too bad,
						 * so go on to the next one
						 */
						++child1It;
						discSol1 = 0;
					} else if(++discSol1 >= maxDiscardedSolutions) {
						GAout << "Warning: The algorithm may be stuck. Try increasing the badSolutionThreshold!" << std::endl;
					}
				} catch(const ::Evaluator::EvaluatorException& ee) {
					if(this->ctrl.verbosity >= VERBOSE) {
						GAout << "Could not evaluate chromosome: " << ee.what() << std::endl;
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
					if((this->evaluator.evaluate(**child2It) > cutoff) || (discSol2 > maxDiscardedSolutions)) {
						if((*child2It)->getFitness() < minFitness) {
							minFitness = (*child2It)->getFitness();
						}

						this->addChromosomeToElite(**child2It);

						/*
						 * The child is no duplicate (or accepted as one) and is not too bad,
						 * so go on to the next one
						 */
						++child2It;
						discSol2 = 0;
					} else if(++discSol2 >= maxDiscardedSolutions) {
						GAout << "Warning: The algorithm may be stuck. Try increasing the badSolutionThreshold!" << std::endl;
					}
				} catch(const ::Evaluator::EvaluatorException& ee) {
					if(this->ctrl.verbosity >= ON) {
						GAout << "Could not evaluate chromosome: " << ee.what() << std::endl;
					}
				}
			}

			if(check_interrupt()) {
				this->interrupted = true;
			}
		}

		/*
		 * Transform the fitness map of the current generation to start at 0
		 * and copy old generation to new generation
		 */		
		sumFitness = this->updateCurrentGeneration(newGeneration, minFitness, false);

		if(this->ctrl.verbosity >= VERBOSE && this->ctrl.verbosity != DEBUG_EVAL) {
			this->printCurrentGeneration();
		}

		discSol1 = 0;
		discSol2 = 0;
	}

	for(ChVecIt it = newGeneration.begin(); it != newGeneration.end(); ++it) {
		delete *it;
	}
}
