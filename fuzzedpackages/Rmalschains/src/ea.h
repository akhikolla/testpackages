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

#ifndef _EA_H

#define _EA_H 1

#include "iea.h"
#include "statistics.h"
#include "cross.h"

using namespace std;


namespace realea {

/**
 * @class EAlgorithm
 * @ingroup realea_common
 *
 * @brief It represents an Evolutionary Algorithm 
 *
 */
class EAlgorithm : public Resetable {
public:
    /**
     *  Constructor.
     *
     *  @param alg algorithm 
     */
    EAlgorithm(IEAlgorithm *alg, ProblemParamPtr problem);

    /**
     * Destructor. Remove the Populations.
     */
    virtual ~EAlgorithm(void);

    /**
     * Set the maximum number of evaluations 
     * 
     * @param maxeval
     */
    void setMaxEval(unsigned int maxeval);

    /**
     * Allow to set the popsize
     *
     * @param popsize size of population
     */
    virtual void setPopsize(unsigned int popsize);

    /**
     * @return the running criterion. 
     */
    RunningPtr getRunning(void);

    /**
     * Apply the algorithm 
     *
     * @param sol best found solution, output
     * @param pfitness fitness of the best solution, output
     *
     * @return the evaluation numbers required for achieve the solution
     */
    unsigned apply(tChromosomeReal &sol, tFitness *pfitness);


    /**
     * Set the class that obtain the statistical information
     *
     * @param stat Statistical class
     */
    void setShow(Statistics *stat);
private:
    /**
     * Set the problem to be solved
     *
     * @param problem the problem a consider
     */
    void setProblem(ProblemParamPtr problem);

    /**
     * Set the evaluation function (it is problem by default)
     *
     * @param eval evaluation class
     */
    void setEval(IEval *eval); 


    /**
     * Set the default size of population, used if setPopsize is not used
     */
    void setDefaultPopsize(void);


protected:
    /**
     * Eval the solution
     *
     * @param newind individual, it is updated with the fitness
     */
//    virtual void eval(tIndividualReal *newind) {
//	m_alg->eval(newind);
//    }
	
protected:
    IEAlgorithm *m_alg; /**< IEAlgorithm */
    Problem *m_problem; 
    IEval *m_eval;
    const unsigned default_popsize; /**< Default Popsize (obtained by setDefaultSize) */
    internal::CrossBinaryPtr m_cross; /**< Crossover operator */
    Statistics *m_stat; /**< Statistical class */
    RunningPtr m_running;
};

typedef EAlgorithm EA;
}

#endif
