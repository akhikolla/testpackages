/**
 * Copyright 2008, Daniel Molina Cabrera <danimolina@gmail.com>
 * 
 * This file is part of software Realea
 * 
 * Realea is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Realea is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _POPULATION_H

#define _POPULATION_H 1

#include "individual.h"
#include "domain.h"
#include "initind.h"
#include "random.h"
#include "define.h"
#include "popobserver.h"
#include <deque>

using namespace std;

namespace realea {

/**
 * @class PopulationReal
 * @ingroup realea_common
 *
 * @brief Stores the individuals (solutions)
 *
 * Stores a group of individuals (its size m_minsize is fixed, but a temporal greater number m_aditionalsize is allowed)
 */
class PopulationReal {
   public:
	/**
         * Reduce the population to the half of the current popsize
         *
         * Aditionally, clear the observers
         */
        void reduceHalf(void);

	/**
 	 * Constructor. 
 	 * @param random generator of the population
 	 * @param max maximum population size.
 	 * @param pob normal population size.
 	 */
	PopulationReal(Random *random, unsigned int max, unsigned int pob);
	virtual ~PopulationReal(void);
	/**
 	 * Reset the population size randomly. 
 	 *
 	 * @param domain domain of individuals.
 	 * @param init the begin position of individuals to be restarted.
 	 * (init = -1 does not keep anyone. init == x keep the individual)
 	 */
	virtual void reset(DomainRealPtr domain, int init=-1);

	/**
  	 * Removes the (m_aditionalsize - m_minsize) worse individuals. 
  	 */
	void removeWorses(void); 

	/**
 	 * Orders randomly the individuals
 	 */
	void random(void); 

	/**
 	 * Sorts the individuals in function of its fitness
 	 */
	void sort(void); 

	/**
	 * @param pos the position
	 * @return the current individual in position pos
	 */
	tIndividualReal *getInd(unsigned int pos);

	/**
	 * Replaces the current individual in position pos with new individual real
	 *
	 * @param pos of the current individual to be removed (> 0 <= size())
	 * @param ind new individual to be inserted
	 */
//	void replace(int pos, const tIndividualReal *ind);
        /**
	 * Replaces the current individual in position pos with new individual real
	 *
	 * @param pos of the current individual to be removed (> 0 <= size())
	 * @param ind new individual to be inserted
	 */
	void replace(unsigned pos, tIndividualRealPtr newind);


	/**
	 * Append a new individual
	 *
	 * @param ind new individual
	 */
	void append(tIndividualReal *ind);
	void append(tChromosomeReal &sol, tFitness fitness);
	/**
 	 *  change the current position with individual sol
 	 */
	void change(unsigned int pos, const tChromosomeReal &sol, tFitness fitness);


	/**
	 * Replaces the current individual in position pos with new individual real
	 * without deleting the current individual in position pos, to allow use it
	 * externally
	 *
	 * @param pos of the current individual to be removed from the
	 * population (> 0 <= size())
	 * @param ind new individual to be inserted
	 */
        void replaceWithoutDeleting(unsigned pos, tIndividualRealPtr newind);

	/**
	 * @return the position of best individual (by fitness)
	 * @see getBests
	 */
	unsigned getBest(void);

	/**
	 * @return the second best fitness from the population
	 */
	tFitness getSecondBestFitness(void);

	/**
	 * @return the Median of the population's fitness
	 */
	tFitness getMedian(void);
	/**
	 * @return the mean of the population's fitness
	 */
	tFitness getMean(void);
	/**
	 * Obtain the difference percentiles (useful to test the behaviour
	 * population 
	 */
	void getPercentils(double *percen, unsigned num); 
	/**
	 * @return the positions of best individuals (by fitness)
	 *
	 * @param num the number of best individuals to return. 
	 * (if num == 1 it is equals that call to getBest)
	 * @see getBest
	 */
	vector<unsigned> getBests(unsigned num);

	/**
	 * @return the position of worst individual (by fitness)
	 */
	unsigned getWorst(void);
	/**
	 * @return the actual size (ignoring individual)
	 */
	unsigned size(void);
	/**
	 * @return dimension of individuals
	 */
	unsigned ndim(void);
	/**
	 * The current new individuals are evaluated.
	 *
	 * @param fitness Fitness class.
	 * @param neweval maximum new evaluations (-1 not limit)
	 */
	void eval(IEvalInd *fitness, unsigned newevals=-1);

	/**
	 * add Observer
	 */
	void setObserver(PopulationObserver *observer);

	/**
         * get a new instance of a chromosome 
         */
	virtual tIndividualReal* getInstance(tChromosomeReal &crom, tFitness fitness);

 	/**
         * get a new instance of a chromosome 
         */
	virtual tIndividualReal* getInstance(tChromosomeReal &crom);

        /**
         * set the population
         */
        void setIndividuals(vector<tIndividualReal*> inds){m_individuals = inds;}


        /*
         * get the diversity of the popualation
         */
        double getDiversity(void);


   private:
	/**
	 * Update the identifications of individuals (it is called after removeWorses)
	 * 
	 * @see removeWorses sort
	 */
	void updateObservers(void); 

	/**
	 * Update the identifications of individuals (it is called after removeWorses)
	 * 
	 * @see removeWorses sort
	 */
	void resetObservers(void);

	/**
	 * @return true if there is an individual null
	 */
	bool thereNull(void);

	/**
	 * Remove all individuals between begin and end (included)
	 *
	 * @param begin initial elem (included)
	 * @param end finish elem (excluded)
	 */
	void remove(unsigned begin, unsigned end); 
	/**
	 * Remove all individuals with not been evaluated
	 */
	void removePending(void);

	unsigned int m_aditionalsize, /**< Aditional size (for storing temporally more individual */
			 m_size; /**< Normal size of population */
	vector<tIndividualReal*> m_individuals; /**< vector of individuals stored into the population */

	unsigned m_worst; /*< position of worst individual */
	unsigned m_best; /*< position of best individual */
	

	bool m_knownbest; 
	bool m_knownworst;
	InitIndividual *m_initInd; /*< Creator of new individuals */
	deque<PopulationObserver *>m_observers; /**< Store the class to notify */

protected:
	/**
	 * Update the identifications of a specific individual
	 *
	 * @param id individual id to notify
	 */
	void notifyObservers(unsigned id); 


	Random *m_random; /*< random generator */
};

typedef PopulationReal* PopulationRealPtr;

}

#endif
