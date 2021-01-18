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

#ifndef _ILOCALSEARCH_H

#define _ILOCALSEARCH_H 1

#include "problem.h"
#include "running.h"
#include "populationreal.h"
#include "random.h"

namespace realea {

/**
 * This class store the internal state of the LS. 
 *
 * @ingroup realea_common
 *
 * each ILocalSearch works with a subclass of that class to store its
 * state (allowing to continue the local search in next running)
 */
class ILSParameters {
public:
virtual ~ILSParameters(void){}
};

/**
 * @class ILocalSearch 
 *
 * @brief This interface allow to define how to obtain the local search method
 *
 * To define new Local Search methods it is only required to implement this interface and
 * define apply and getInitOptions methods.
 */
class ILocalSearch {
public:
    /**
     * Obtain the initial parameter to local searchObtain the initial parameter to local search.
     *
     * @param sol chromosome to be improved
     * @see applyNewSol
     */
    virtual ILSParameters *getInitOptions(tChromosomeReal &sol)=0;

    /**
      *
      * Apply the LS method improvement 
      *
      * @param params the initial parameters of the LS, initialised by getInitOptions
      * @param sol chromosome to be improved, it is changed
      * @param fitness fitness of chromosome 'sol', it is updated by the new sol values
      * @param nstep evaluation number into the LS process
      * @return evaluations number
      *
      * @see getInitOptions
      */
     virtual unsigned apply(ILSParameters *params, tChromosomeReal &sol, tFitness &fitness, unsigned nstep)=0;

     /**
      * Allow to specify the population whose individuals will be improved
      *
      * @param pop population
      */
     void setPopulation(PopulationReal *pop) {
	m_pop = pop;
     }

     /**
      * Set the running criterion
      *
      * @param run stop criterion
      */
     void setRunning(RunningPtr run) {
	m_running = run;
     }

     /**
      * Set the problem 
      *
      * @param problem problem to be evaluated
      */
     void setProblem(Problem *problem) {
	m_problem = problem;
     }

     /**
      * Set the Random number generator
      *
      * @param random 
      */
     void setRandom(Random *random) {
	m_random = random;
     }

     void setEval(IEval *eval) {
	m_eval = eval;
     }

     virtual ~ILocalSearch(void) {}

    
    protected:
     Random *m_random; /**< The current randomgeneration numbers */
     IEval *m_eval; /** The current evaluation function */
     Problem *m_problem; /**< The current problem */
     PopulationReal *m_pop; /** The current population */
     RunningPtr m_running; /**< The stopping criterion */
};

/**
 * This class allow us to recover and receive the parameters using an array (for Niching Techniques)
 */
class IParallelLocalSearch : public ILocalSearch {
public:
    /**
     * Obtain the parameter from a serialized array.
     *
     * @param params array of params
     * @param size size of array
     */
    virtual ILSParameters *recoverOptions(tGen *params, unsigned size)=0;

    /**
     * Store in a real vector the values of a parameters. 
     *
     * @param params parameter to store
     * @param paparams reference to sequencial parameters (can be NULL)
     * @param size reference to size (always is right, even if params == NULL)
     *
     * @see recoverOptions
     */
    virtual void storeOptions(ILSParameters *params, tGen **paparams, unsigned *psize)=0;
};


}

#endif
