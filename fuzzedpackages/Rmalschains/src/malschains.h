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
#ifndef _MALSChains_H

#define _MALSChains_H 1

#include "restart.h"
#include "selectls.h"
#include "reducepopulation.h"
#include "lsparammem.h"

#include "hybrid.h"

namespace realea {

/**
 * @class MALSChains
 *
 * @brief memetic algorithm following the model proposed in the Daniel Molina's thesis
 */
class MALSChains : public Hybrid {
public:
	/**
	 * Constructor
	 *
	 * 
	 */
	MALSChains(IEAlgorithm *alg, ILocalSearch *ls);

	~MALSChains(void);

	/**
	 * Set the Threshold for restart the population. When a certain times neither the EA and the LS
	 * can improve the solution, the population is restarted (keeping the best value). 
	 *
	 * @param threshold minimum ratio between improvement and fitness. If the minimum is lower than 
	 * this value, a failed is register.
	 *
	 * @param maxfailed maximum of fails to restart
	 */
	void setThreshold(double threshold, int maxfailed);
	/**
	 * @param ratio. Set the global ratio invested into the Local Search
	 *
	 * @param ratio global ls/total ratio
	 * (change the intensity, because only one LS improvement is made)
	 */
	void setEffortRatio(double ratio);
	void setDebug(void);

	void recoverIndividual(unsigned oldind, tGen *aind, unsigned size, tGen *aoptional, unsigned size_optional);
	void storeIndividual(tIndividualRealPtr ind, tGen **paind, unsigned *pmax, tGen **padditional, unsigned *pmaxad); 

	/**
	 * Enable the population size reduction, using PopulationReductionStrategy
	 *
	 * @param reduction number
	 */
	void enablePopReduction(unsigned numReduction);
	void setMaxEval(unsigned int maxeval);

	void setSelectImprovementStrategy(SelectImprovementLS *select_improvement) {
	    m_select_improvement = select_improvement;
	}

	void setRunning(Running *running); 

	RunningPtr getRunning(); 

        unsigned getNumEvalEA();
        unsigned getNumEvalLS();

	int getNumImprovementEA();
	int getNumImprovementLS();
	int getNumTotalEA();
	int getNumTotalLS();
	tFitness getImprovementEA();
	tFitness getImprovementLS();
	double getTimeMsEA();
	double getTimeMsLS();
	double getTimeMsMA();

	/**
	 * Set the restart strategy. If it is not called never it will make a restart.
	 */
	void setRestart(RestartStrategy *restart);

	/**
	 * Alter the EA and after the LS, then the EA again, the LS, ...
	 */
	unsigned realApply(tChromosomeReal &sol, tFitness &fitness);
	unsigned init(void);

	void setInitEval(IEval* eval); 

	/**
	 * Set a distorsion size when the solution is a local optimum
	 *
	 * @param size size of maximum disruption
	 */
	void setDisruptionSize(double size); 
protected:
	/**
	 * Check if there is enough diversity
	 */
	bool hasDiversity(PopulationReal *pop);

	/**
	 * Set the differences
	 */
	void setDif(bool debug, string ident, unsigned id, tFitness oldfit, tFitness newfit);

protected:
	/**
	 * Disturb the solution (similar to mutation
	 * it uses the disruption_size
	 *
	 */
	void disturb(tChromosomeReal &sol);

	/**
	 * Check the improvemt was relevant
	 */
	bool hasImprovedEnough(tFitness oldfitness, tFitness fitness); 
private:
	PopulationReductionStrategy m_popReductions;

protected:
	double m_disruption_size;
	tFitness m_threshold;
	internal::LSParametersMemory *m_memory;
	unsigned m_nevalalg;
	unsigned m_nevalls;
	unsigned m_maxfailed;
	IEvalInd *m_initeval;
	unsigned m_initMaxEval;
	double m_effort;
	RestartStrategy *m_restart;
	SelectImprovementLS *m_select_improvement;
	bool m_debug;

    int m_num_improvement_ea, m_num_improvement_ls;
    int m_num_total_ea, m_num_total_ls;
    tFitness m_improvement_alg, m_improvement_ls;
    double m_time_ms_alg, m_time_ms_ls, m_time_ms_ma;


};


}

#endif
