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

#include "malschains.h"
#include "selectlsimp.h"
#include "get_util.h"
#include "debug.h"
#include <cassert>
#include <cstdio>
#include <ctime>
#include <cmath>

#define _DEBUG_POP 0

using namespace realea;
using namespace realea::internal;

void MALSChains::setEffortRatio(double ratio) {
	if (ratio == 1)
	   throw new string("MALSChains::effortRatio is not valide");

	m_effort = ratio;
}

void MALSChains::setThreshold(double threshold, int maxfailed) {
    m_threshold = threshold;
    m_maxfailed = maxfailed;
    
    print_info("Threshold: %d[%Le]\n", m_maxfailed, m_threshold);
}

MALSChains::MALSChains(IEAlgorithm *alg, ILocalSearch *ls) : Hybrid(alg,ls), m_memory(NULL) {
   m_effort = -1;
   m_nevalalg = m_nevalls = 0;
   m_restart = NULL;
   m_debug = false;
   m_select_improvement = NULL;
   m_disruption_size = 0;
   m_threshold = 0;
   m_maxfailed = 0;
}

void MALSChains::setDisruptionSize(double size) { 
   assert(size >= 0 && size < 1);
   m_disruption_size = size;
}


void MALSChains::setDebug(void) {
    m_debug = true;
    enable_print_debug();
}

MALSChains::~MALSChains(void) {
    if (m_memory) {
	delete m_memory;
    }

    if (m_restart) {
	delete m_restart;
    }

    if (m_select_improvement) {
	delete m_select_improvement;
    }

    if (m_initeval) {
	delete m_initeval;
    }
}

/**
 * Calculate the new frec to apply to obtain the indicated ratio.
 *
 * It is required because there are AEs that could not follow the frec indicated in the previous step
 *
 * @param nevalalg total number evaluation during the EA
 * @param nevalls total number evaluation during the LS
 * @param intensity LS intensity
 *
 * @return frec to use the EA to maintain the ratio
 */
unsigned calculateFrec(unsigned nevalalg, unsigned nevalls, unsigned intensity, double ratio) {

    assert(ratio > 0);

    double coc = ratio*(nevalalg + nevalls + intensity) - nevalalg;
    double div = (1-ratio);

    if (div == 0) {
	throw string("MALSChains::Ratio is too low");
    }

    return ((unsigned) floor(coc/div));
}

void MALSChains::enablePopReduction(unsigned numReduction) {
    m_popReductions.setNumReductions(numReduction);
}

unsigned MALSChains::init(void) {
    initLs();

    m_initMaxEval = m_running->maxEval();
    m_popReductions.config(m_running); 
    unsigned neval = m_alg->init();

    if (m_select_improvement == NULL) {
	m_select_improvement = new SelectBestToImprove();
    }

    if (m_memory == NULL) {
	m_memory = new LSParametersMemory(m_alg->getPop()->size());
	m_alg->getPop()->setObserver(m_memory);
    }

    m_nevalalg = neval;
    m_nevalls = 0;
    return neval;
}

void MALSChains::setMaxEval(unsigned int maxeval) {
    assert(maxeval >= m_intensity);
    unsigned frec = calculateFrec(m_nevalalg, m_nevalls, m_intensity, m_effort);
    m_alg->setMaxEval(frec);
}

void MALSChains::setRestart(RestartStrategy *restart) {

    if (m_restart) {
	delete m_restart;

    }

    m_restart = restart;
}

RunningPtr MALSChains::getRunning(void) {
    return m_alg->getRunning();
}

unsigned MALSChains::getNumEvalEA(void) {
    return m_nevalalg;
}

unsigned MALSChains::getNumEvalLS(void) {
    return m_nevalls;
}

int MALSChains::getNumImprovementEA(void) {
    return m_num_improvement_ea;
}

int MALSChains::getNumImprovementLS(void) {
    return m_num_improvement_ls;
}

int MALSChains::getNumTotalEA(void) {
    return m_num_total_ea;
}

int MALSChains::getNumTotalLS(void) {
    return m_num_total_ls;
}

tFitness MALSChains::getImprovementEA(void) {
    return m_improvement_alg;
}

tFitness MALSChains::getImprovementLS(void) {
    return m_improvement_ls;
}

double MALSChains::getTimeMsEA(void) {
    return m_time_ms_alg;
}

double MALSChains::getTimeMsLS(void) {
    return m_time_ms_ls;
}

double MALSChains::getTimeMsMA(void) {
    return m_time_ms_ma;
}


void MALSChains::setRunning(Running *running) {
    ProxyEA::setRunning(running);
    unsigned frec = calculateFrec(m_nevalalg, m_nevalls, m_intensity, m_effort);
    m_alg->setRunning(m_running->getSubRunning(frec));
}



void MALSChains::recoverIndividual(unsigned oldind, tGen *aind, unsigned size, tGen *aoptional, unsigned size_optional) {
    m_alg->recoverIndividual(oldind, aind, size-1, aoptional, size_optional);

    if (aind[size]) {
        m_alg->getPop()->getInd(oldind)->incremCount("non_improved");
    }

    if (aoptional != NULL) {
        IParallelLocalSearch *ls = (IParallelLocalSearch *) m_ls;
	ILSParameters *params = ls->recoverOptions(aoptional, size_optional);
	assert(m_memory);
	m_memory->store(oldind, params);
    }

}

void MALSChains::storeIndividual(tIndividualRealPtr ind, tGen **paind, unsigned *pmax, tGen **paoptional, unsigned *psize_optional) {
    tGen *asol, *aoptional_sol;
    tGen *asol_ma;
    unsigned size_sol, size_optional;
    tGen *aparams;
    unsigned size_param;
    unsigned size;

    m_alg->storeIndividual(ind, &asol, &size_sol, &aoptional_sol, &size_optional);
    assert(size_optional == 0 && aoptional_sol == NULL);

    asol_ma = new tGen[size_sol+1];
    copy(asol, asol+size_sol, asol_ma);
    delete[] asol;
    asol_ma[size_sol] = (ind->getCount("non_improved") > 0) ? 1 : 0;

    *paind = asol_ma;
    *pmax = size_sol+1;

    size_param = 0;
    
    if (m_memory) {
        unsigned posind = ind->getId();
        IParallelLocalSearch *ls = (IParallelLocalSearch *) m_ls;
	ILSParameters *params = m_memory->recover(posind);
	ls->storeOptions(params, &aparams, &size_param);
    }

    size = size_optional+size_param;
    assert(size > 0);

    *psize_optional = size_param;
    *paoptional = NULL;

    if (aoptional_sol != NULL || aparams != NULL) {
       *paoptional = new tGen[size];

       if (aoptional_sol != NULL) {
	  copy(aoptional_sol, aoptional_sol+size_optional, *paoptional);
	  delete[] aoptional_sol;
       }

       if (aparams != NULL) {
	  copy(aparams, aparams+size_param, *paoptional+size_optional);
	  delete[] aparams;
       }

    }
}

void MALSChains::setDif(bool debug, string ident, unsigned id, tFitness oldfit, tFitness newfit) {
   if (debug) {

      if (oldfit!= newfit) {
     print_info("%s[%2d]:\t%e -> %e  %e\n", ident.c_str(), id, oldfit, newfit, fabs(newfit-oldfit));
      }
//      else {
//	 print_info("%s[%2d]:\t%Le\n", ident.c_str(), id, newfit);
//      }
   }

}

bool MALSChains::hasDiversity(PopulationReal *pop) {
    return true;
    double percen[5];
    pop->getPercentils(percen, 4);


	    print_info("EA::Improvement: %e\t%e\t%e\t%e\t%e\n", percen[0], percen[1], 
			percen[2], percen[3], percen[4]);


    if (percen[2] == percen[4]) {
	return false;
    }
    else if (percen[1] == percen[3]) {
	return false;
    }
    else if (fabs((percen[0] - percen[2])/percen[2]) < 1e-3) {
	return false;
    }
    else {
	return true;    
    }
}

void MALSChains::setInitEval(IEval*eval) {
   Hybrid::setInitEval(eval);
   m_initeval = (IEvalInd *) new EvalRunning(eval, m_running);
}

void MALSChains::disturb(tChromosomeReal &sol) {
    DomainRealPtr domain = m_problem->getDomain();
    unsigned dim = domain->getDimension();
    double min, max;
 
    for (unsigned i = 0; i < dim; i++) {
        if (domain->canBeChanged(i)) {
	    domain->getValues(i, &min, &max);
	    sol[i] += m_disruption_size*m_random->randreal(-1,1)*(max-min);
	}
    }

    domain->clip(sol);
}

bool MALSChains::hasImprovedEnough(tFitness oldfitness, tFitness fitness) {
   bool hasImproved;
   double minimum = m_running->getThreshold()/10.0;

   // Check the solution has improvement and if m_threshold is defined, the fitness/improvement is
   // greater that value
   if (!m_problem->isBetter(fitness, oldfitness)) {
      hasImproved = false;
   }
   else if (fabs(fitness-oldfitness) < minimum) {
      hasImproved = false;
   }
   else if (m_threshold == 0) {
      hasImproved = true;
   }
   else {
      hasImproved = (fabs(oldfitness-fitness)/fabs(fitness) >= m_threshold);

//      if (m_debug && !hasImproved) {
//	 print_info("Not Improvement: %Le -> %Le is lower than %Le: %Le\n", oldfitness, fitness, m_threshold, fabs(oldfitness-fitness)/fabs(fitness));
//      }
   }

   return hasImproved;
}

void printPopFitness(double *fitnessold, double *fitness, unsigned num) {
//    for (unsigned i=0; i < num; i++) {
//	print_info("PopFitness[%d]: %e->%e\n", i, fitnessold[i], fitness[i]);
//    }

#ifdef _DEBUG_POP   
    print_info("EA::PopFitness:  ");

    for (unsigned i=0; i < num; i++) {
	print_info(" %e ", fitness[i]);
    }

    print_info("\n");


    print_info("EA::Improvement: ");

    for (unsigned i=0; i < num; i++) {
	print_info(" %e ", fabs(fitnessold[i]-fitness[i]));
    }

    print_info("\n");
#endif
}

unsigned MALSChains::realApply(tChromosomeReal &bestsol, tFitness &bestfitness) {
    tIndividualReal *ind, *best;
    unsigned posind;
    tFitness oldfitness, fitness, fitness_alg;
    unsigned ndim = bestsol.size();
    tChromosomeReal sol(ndim), sol_alg(ndim);
    unsigned alg_failed;
    clock_t m_time_ls, m_time_alg, m_time_ma;
    clock_t clock_begin, m_time_ma_begin;
    bool hasImprovedByEA,hasImprovedByLS;
    unsigned restarts=0;
    
    PopulationReal *pop_alg = m_alg->getPop();
    deque<tIndividualReal*> ind_to_improve;

    m_time_ls = m_time_alg = m_time_ma = 0;
    m_num_improvement_ea = m_num_improvement_ls = 0;
    m_num_total_ea = m_num_total_ls = 0;

    sol=bestsol;
    fitness = bestfitness;
    m_improvement_alg = m_improvement_ls = 0;

    unsigned initMax = m_running->numEval();
    fitness_alg = pop_alg->getInd(pop_alg->getBest())->perf();

    alg_failed = 0;
    m_time_ma_begin = clock();

//    while ( (m_nevalalg+m_nevalls) < m_initMaxEval && !m_running->hasFoundOptime()) {
    while (!m_running->isFinish()) {
        double medidas[5], medidasold[5];
        tFitness old_secondfitness = pop_alg->getSecondBestFitness();
        tFitness old_fitness = fitness_alg;
	clock_begin = clock();

	if (m_debug) {
	    pop_alg->getPercentils(medidasold, 4);
	}

	m_nevalalg += m_alg->realApply(sol_alg, fitness_alg);

	if (m_debug) {
	    pop_alg->getPercentils(medidas, 4);
	    printPopFitness(medidasold, medidas, 4);
	}

	m_time_alg += clock()-clock_begin;
	m_improvement_alg += fabs(fitness_alg-old_fitness);
        tFitness new_secondfitness = pop_alg->getSecondBestFitness();
	hasImprovedByEA = hasImprovedEnough(old_secondfitness, new_secondfitness);
	assert(new_secondfitness <= old_secondfitness);

	if (m_debug) {
	   if (restarts) {
	      setDif(m_debug, "EASBest", -1, new_secondfitness, old_secondfitness);
	   }
	   else {
	      setDif(m_debug, "EABest", pop_alg->getBest(), old_fitness, fitness_alg);
	//      setDif(m_debug, "EASecondBest", -1, old_secondfitness, new_secondfitness);
	   }
	}

	if (fitness_alg != old_fitness) {
	    m_num_improvement_ea++;
	}
	m_num_total_ea++;

	// Check the optime
	if (m_running->isOptime(fitness_alg)) {
	    continue;
	}

	// Select the new individual to increm
	m_select_improvement->getIndsToImprove(pop_alg, ind_to_improve);

	if (ind_to_improve.size()!=0) {
	    // Select the individual to improve
	    posind = m_select_improvement->selectIndToImprove(ind_to_improve);
	}
	else {
	    // Choose the new one randomly
	    posind = m_random->randint(0, pop_alg->size()-1);
	}

	ind = pop_alg->getInd(posind);
	sol = ind->sol();

	ILSParameters *params = m_memory->recover(posind);
	bool recover = (params != NULL);
	
	if (params == NULL) {	
	    // Apply the LS to the best one with the rest of intensity
	    params = m_ls->getInitOptions(sol);
	}
	
	fitness = ind->perf();
	oldfitness = fitness;

	clock_begin= clock();
        m_nevalls += 
	    m_ls->apply(params, sol, fitness, m_intensity);
	m_time_ls += clock()-clock_begin;

	ind->incremCount("ls"); 
	m_improvement_ls += fabs(fitness-oldfitness);

	setDif(m_debug, "LS ", ind->getId(), oldfitness, fitness);

	hasImprovedByLS = hasImprovedEnough(fitness_alg, fitness);

	if (m_problem->isBetter(fitness, fitness_alg)) {
	    fitness_alg = fitness;
	    m_num_improvement_ls++;
	}
	
	
	if (!hasImprovedEnough(oldfitness, fitness)) {
	   ind->incremCount("non_improved");

	   if (recover) {
	      m_memory->remove(posind);
	   }
	   else {
	      delete params;
	   }

//         TODO: fix the disruption_size operator
//	   if (m_disruption_size > 0) {
//	      unsigned pos_best = pop_alg->getBest();
//	      tIndividualRealPtr ind_dis=NULL;
//
//	      if (posind != pos_best) {
//		 ind_dis = ind;
//	      }
//	      else {
//		 if (lastpos_better >= 0) {
//		 ind_dis = pop_alg->getInd(lastpos_better);
//
//		 if (ind_dis->perf() != lastpos_fitness || pos_best == lastpos_better) {
//		    ind_dis = NULL;
//		 }
//		}
//
//		 lastpos_better = pos_best;
//		 lastpos_fitness = pop_alg->getInd(lastpos_better)->perf();
//	      }
//
//	      if (ind_dis != NULL) {
//		 tChromosomeReal sol = ind_dis->sol();
//		 disturb(sol);
//		 fitness = m_eval->eval(sol);
//		 pop_alg->change(ind_dis->getId(), sol, fitness);
//	      }
//
//	   }

	   
//	   if (m_disruption_size == 0 && ind_to_improve.empty() && m_restart != NULL && !m_running->isFinish()) 
	}
	else {
	    pop_alg->change(posind, sol, fitness);
	    m_memory->store(posind, params);
	}

	if (!hasImprovedByEA && !hasImprovedByLS) {
	   alg_failed++;
	}
	else {
	   alg_failed = 0;
	}

	m_num_total_ls++;

	bool restart = false;

	if (ind_to_improve.size() == 0 && !hasImprovedByLS) {
	   restart = true;
	}

	if ((m_maxfailed != 0 && alg_failed >= m_maxfailed && !m_running->isFinish())) 
	   restart = true;

	if (m_restart != NULL && restart) {
		m_restart->apply(pop_alg, m_problem, m_initeval);
		restarts++;
		alg_failed = 0;
		
		print_info("Restart_AlgFailed\t: %Le\n", fitness_alg);

		fitness_alg = pop_alg->getInd(pop_alg->getBest())->perf();

		print_info("Restart_AlgNewBest\t: %Le\n", fitness_alg);

		m_memory->reset();
	}

	// Update the population size
	if (m_popReductions.updatePopulationSize(pop_alg)) {
	   print_info("ReducedPopulation: %u\n", m_running->numEval());
	}

    }

    m_time_ma = clock()-m_time_ma_begin;

    if (m_debug) {
        double ratio_effort = ( (double) m_nevalalg)/(m_nevalalg+m_nevalls);
	print_debug("RatioEffort Alg/LS: [%.0f/%.0f]\n", 100*ratio_effort, 100*(1-ratio_effort));

        double ratio_alg = m_improvement_alg/(m_improvement_alg+m_improvement_ls);
	print_debug("RatioImprovement Alg/LS: [%.0f/%.0f]\n", 100*ratio_alg, 100*(1-ratio_alg));
	print_debug("Restarts: %d\n", restarts);
    }
    


    m_time_ms_alg = (m_time_alg*1000.0)/CLOCKS_PER_SEC;
    m_time_ms_ls = (m_time_ls*1000.0)/CLOCKS_PER_SEC;
    m_time_ms_ma = (m_time_ma*1000.0)/CLOCKS_PER_SEC;

    print_debug("Time[ALG]: %.2f\n", m_time_ms_alg);
    print_debug("Time[LS]: %.2f\n", m_time_ms_ls);
    print_debug("Time[MA]: %.2f\n", m_time_ms_ma);
    print_debug("RatioTime[ALG/MA]: %.2f\n", 100*m_time_ms_alg/m_time_ms_ma);
    print_debug("RatioTime[LS/MA]: %.2f\n", 100*m_time_ms_ls/m_time_ms_ma);
 //   print_debug("RatioTime[(ALG+LS)/MA]: %.2f\n", 100*(m_time_ms_alg+m_time_ms_ls)/m_time_ms_ma);

    if (m_num_total_ea > 0) {
       print_debug("NumImprovement[EA]:%d%%\n", (m_num_improvement_ea*100)/m_num_total_ea);
    }
    if (m_num_total_ls > 0) {
       print_debug("NumImprovement[LS]:%d%%\n", (m_num_improvement_ls*100)/m_num_total_ls);
    }
    
    unsigned neval = m_running->numEval()-initMax;
    best = pop_alg->getInd(pop_alg->getBest());
    bestsol = best->sol();
    bestfitness = best->perf();
    m_running->reset();
    return neval;
}
