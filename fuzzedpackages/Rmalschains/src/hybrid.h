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

#ifndef _HYBRID_H

#define _HYBRID_H 1

#include "iea.h"
#include "random.h"
#include "ilocalsearch.h"

namespace realea {

/**
 * @class ProxyEA
 *
 * @brief Allow to proxy an EAlgorithm (using the pattern decorator)
 *
 * This class contains an EA, and redirect all method to this variable
 */
class ProxyEA : public IEAlgorithm {
public:
     ProxyEA(IEAlgorithm *alg) : m_alg(alg) {}

     unsigned getDefaultPopsize(void) {
	 return m_alg->getDefaultPopsize();
     }

     virtual unsigned realApply(tChromosomeReal &sol, tFitness &fitness)=0;
     virtual unsigned init(void)=0;

    /**
     * Allow to set the popsize
     *
     * @param popsize size of population
     */
    void setPopsize(unsigned int popsize) {
	m_alg->setPopsize(popsize);
    }

    void setProblem(Problem *problem) {
	m_problem = problem;
	m_alg->setProblem(problem);
    }

    void setRunning(Running *running) {
	m_running = running;
    }

    virtual void recoverIndividual(unsigned pos, tGen *aind, unsigned size, tGen *aoptional, unsigned addsize) {
	m_alg->recoverIndividual(pos,aind, size, aoptional, addsize);
    }

    virtual void storeIndividual(tIndividualRealPtr ind, tGen **paind, unsigned *pmax, tGen **paoptional, unsigned *paddsize) {
	m_alg->storeIndividual(ind, paind, pmax, paoptional, paddsize);
    }

    virtual void setRandom(Random *random) {
	m_alg->setRandom(random);
    }

    virtual Random *getRandom(void) {
	return m_alg->getRandom();
    }

    void setStat(Statistics *stat) {
	m_alg->setStat(stat);
    }

    virtual void setMaxEval(unsigned int maxeval)= 0;
 
    unsigned getMaxEval(void) {
	return m_alg->getMaxEval();
    }

    void reset(void) {
    	m_alg->reset();
    }

    PopulationReal *getPop(void) {
	return m_alg->getPop();
    }	

    virtual ~ProxyEA(void) {
	if (m_alg)
	   delete m_alg;

	if (m_eval)
	    delete m_eval;
    }

    	virtual void setInitEval(IEval*eval) {
		m_alg->setInitEval(eval);
    	}

    	virtual void setNewEval(IEval*eval); 


protected:
	IEAlgorithm *m_alg;
	Problem *m_problem;
	IEvalInd *m_eval;
	Running *m_running;
};

/**
 * @class Hybrid
 *
 * @brief Hybrid scheme, it combines an EA and a LS method (with a Intensity fixed).
 */
class Hybrid : public ProxyEA {
	public:
		/**
		 * Constructor
		 *
		 * @param alg EA to applied
		 * @param ls Local Search to apply
		 */
		Hybrid(IEAlgorithm *alg, ILocalSearch *ls);
 
		/**
		 * @param ratio. Set the global ratio invested into the Local Search
		 *
		 * @param ratio global ls/total ratio
		 *
		 */
		virtual void setEffortRatio(double ratio)=0;

		/**
		 * Set the intensity
		 * @param intensity
		 */
		void setIntensity(unsigned intensity) {
		     m_intensity = intensity;
		}

		~Hybrid(void) {
			if (m_ls) {	
				delete m_ls;	
			}
		}

	/**
	 * init the LS method
	 */
	void initLs(void); 

    virtual void setRandom(Random *random) {
	m_alg->setRandom(random);
	m_random = random;
    } 

    Random *getRandom(void) {
	return m_alg->getRandom();
    } 



protected:
	ILocalSearch *m_ls; /**< LS Method */
	unsigned m_intensity;  /**< LS intensity to apply in each improvement */
	Random *m_random;	
};


}

#endif
