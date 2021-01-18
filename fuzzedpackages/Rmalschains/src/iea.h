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

#ifndef _IEA_H

#define _IEA_H 1

#include "problem.h"
#include "random.h"
#include "populationreal.h"
#include "statistics.h"
#include "cross.h"

using namespace std;

namespace realea {

namespace internal {
class EvalRunning : public IEvalInd {
public:
   EvalRunning(IEval *eval, Running *m_running): m_eval(eval), m_running(m_running) {
   }

   tFitness eval(const tChromosomeReal &sol) {
	tFitness fit = m_eval->eval(sol);
	m_running->notifyEval(fit);
	return fit;
   }

   unsigned eval(tIndividualReal *newind) {
	newind->eval(this);
	return 1;
   }
   ~EvalRunning(void) {}

private:
       IEval *m_eval;  
       Running *m_running;
};
}


class IEAlgorithm  : public Resetable {
public:
     /**
      * @return the default popsize
      */
    virtual unsigned getDefaultPopsize(void)=0;

    /**
     * Init the algorithm before to call realApply
     */
    virtual unsigned init(void)=0;

	/**
	 * Real implementation that applies the algorithm 
	 *
	 * @param sol best found solution, output
	 * @param pfitness fitness of the best solution, output
	 */
	virtual unsigned realApply(tChromosomeReal &sol, tFitness &fitness)=0;

    /**
     * Allow to set the popsize
     *
     * @param popsize size of population
     */
    virtual void setPopsize(unsigned int popsize)=0;

    virtual void setRandom(Random *random)=0;

    virtual Random *getRandom(void)=0;

    virtual void setProblem(Problem *problem)=0;

    virtual void setRunning(Running *running)=0;
    virtual RunningPtr getRunning(void)=0;
    virtual void setStat(Statistics *stat)=0;

    /**
     * Set the maximum eval
     *
     * @param maxeval maximum evaluations number
     */
    virtual void setMaxEval(unsigned int maxeval)=0; 
    virtual unsigned getMaxEval(void)=0;

    /**
     * Eval the solution
     *
     * @param newind individual, it is updated with the fitness
     */
    virtual void reset(void)=0;

    virtual PopulationReal *getPop(void)=0; 	

    virtual void setNewEval(IEval*eval)=0; 

    virtual void setInitEval(IEval*eval)=0; 


/**
 * Stores an individual in a vector (to store in a file or to send to other computers)
 *
 * @param ind Individual to be stored
 * @param paind reference to array to be stored (it is created), it is fixed
 * @param pmax reference to size of the array
 * @param paditional reference to array of aditional information (it can be NULL)
 * @param pmaxad reference to size of the aditional array (0 if there is none)
 *
 * @see recoverIndividual
 */
virtual void storeIndividual(tIndividualRealPtr ind, tGen **paind, unsigned *pmax, tGen **paditional, unsigned *pmaxad)=0;

/**
 * Load an individual from a sequence representation in the position indicated 
 *
 * @param posind position to be replaced
 * @param aind reference to array when it is stored
 * @param reference to size size of the array
 *
 * @return individual recover
 *
 * @see storeIndividual
 */
virtual void recoverIndividual(unsigned posind, tGen *aind, unsigned size, tGen *aoptional, unsigned size_optional) = 0;
};

typedef IEAlgorithm IEA;


class ClassEAlgorithm  : public IEAlgorithm {
public:


     /**
      * @return the default popsize
      */
    virtual unsigned getDefaultPopsize(void)=0;

	/**
	 * Real implementation that applies the algorithm 
	 *
	 * If it is applied several times it must continue from last time
	 *
	 * @param sol best found solution, output
	 * @param pfitness fitness of the best solution, output
	 *
	 */
	virtual unsigned realApply(tChromosomeReal &sol, tFitness &fitness)=0;

    /**
     * Allow to set the popsize
     *
     * @param popsize size of population
     */
    virtual void setPopsize(unsigned int popsize);

	ClassEAlgorithm(Random *random);

	void setRandom(Random *random) {
	     m_random = random;
	}

	Random *getRandom(void) {
	     return m_random;
	}


	void setProblem(Problem *problem) {
	    m_problem = problem;
	}

	void setRunning(Running *running); 
	RunningPtr getRunning(void);

    /**
     * Set the class that obtain the statistical information
     *
     * @param stat Statistical class
     */
    void setStat(Statistics *stat) {
	m_stat = stat;	
    }

    /**
     * Set the maximum eval
     *
     * @param maxeval maximum evaluations number
     */
    void setMaxEval(unsigned int maxeval); 


    unsigned getMaxEval(void); 


    virtual ~ClassEAlgorithm(void);
	
    virtual void setNewEval(IEval*eval); 

    virtual void setInitEval(IEval*eval); 

    virtual void reset(void);

    virtual PopulationReal *getPop(void) {
	 return m_pop;
    }

/**
 * Stores an individual in a vector (to store in a file or to send to other computers)
 *
 * @param ind Individual to be stored
 * @param paind reference to array to be stored (it is created), it is fixed
 * @param pmax reference to size of the array
 * @param paditional reference to array of aditional information (it can be NULL)
 * @param pmaxad reference to size of the aditional array (0 if there is none)
 *
 * @see recoverIndividual
 */
virtual void storeIndividual(tIndividualRealPtr ind, tGen **paind, unsigned *pmax, tGen **paoptional, unsigned *psize_optional);

virtual void recoverIndividual(unsigned pos, tGen *aind, unsigned size, tGen *aoptional, unsigned size_optional);

protected:
	Problem *m_problem; /**< Problem to abord */
	IEvalInd *m_init_eval; /** Eval function */
	IEvalInd *m_new_eval; /** Eval function */

    	RunningPtr m_running; /**< Running criterion */
    	PopulationReal *m_pop; /**< Main Population */
	Random *m_random; /**< Random generator */
	Statistics *m_stat; /**< Statistical class */
	unsigned m_maxeval; 
};


class ICrossEAlgorithm  : public ClassEAlgorithm {
public:
    ICrossEAlgorithm(Random *random); 

    /**
     * Set the crossover to be applied
     *
     * @param cross the binary crossover
     */
    virtual void setCross(ICrossBinaryPtr cross);

    virtual ~ICrossEAlgorithm(void); 

    virtual void setProblem(Problem *problem); 
    virtual void reset(void);

protected:
    internal::CrossBinaryPtr m_cross; /**< Crossover operator */
    ICrossBinaryPtr m_icross;
};

}
#endif
