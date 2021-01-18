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

#include "jdemc.h"
#include <cstdio>
#include <cstring>
#include "random.h"
#include "populationreal.h"
#include "debug.h"
#include <cmath>
#include <cassert>
#include <cstdlib>

#define CopyVector(a,b) memcpy((a),(b),nDim*sizeof(double))
#define Element(a,b,c)  (a->getInd(b).gen(c))

using namespace realea;


JDEMC::JDEMC(Random *random, int reduced) : ClassEAlgorithm(random) {
    NumReduce = reduced;    // 7.2.2010 pmax := NumReduce+1 ==> pmax = 3+1 = 4
    strategy_prob = new double[JDE_NB_OF_STRATEGIES];
    m_strategies[0] = "jDEbin";
    m_strategies[1]= "jDEexp";
    m_strategies[2] = "jDEbest";
    m_debug = false;
}

JDEMC::~JDEMC(void) {
    delete[] strategy_prob;
}

void JDEMC::setDebug(void) {
    m_debug = true;
}

unsigned JDEMC::init(void) {
    // init/check the parameters
    int dim = m_problem->getDimension();
    extDim = dim + 30; // extended dim to store control paramaters

    // Init the population
    m_pop->reset(m_problem->getDomain());
    // Init the different elements: Running, Cross ...
    reset();
    // Eval the population
    m_pop->eval(m_init_eval);
    // Init the counters
    init_strategies = false;
    memset(failure_memory, 0, sizeof(int)*JDE_NB_OF_STRATEGIES*JDE_LP);
    memset(success_memory, 0, sizeof(int)*JDE_NB_OF_STRATEGIES*JDE_LP);
    return m_running->numEval();
}


void printStrategiesStat(string names[], double *strategy_prob) {
    print_info("Probabilities: ");

    for (int k=0; k < JDE_NB_OF_STRATEGIES; k++) {
	print_info("%s[%.2f]", names[k].c_str(), strategy_prob[k]);
    }
    print_info("\n");
}

unsigned JDEMC::realApply(tChromosomeReal &sol, tFitness &fitness) {
    const int D=m_problem->getDimension();
    tChromosomeReal crom(D);
    tChromosomeReal best(D);
    tIndividualRealJDEMC* nov;
    tChromosomeReal bestsol;
    tFitness best_fit;
    unsigned best_ind;
    unsigned it, i;
    tIndividualReal* ind;

    unsigned MaxFES = m_running->maxEval();
    m_running->reset();
    int NP = m_pop->size();
    best_ind = m_pop->getBest();
    best_fit = m_pop->getInd(best_ind)->perf();
    best = m_pop->getInd(best_ind)->sol();


   // While is not finish
   for (it = 0; !m_running->isFinish(); it++) {
        if (m_stat)
	    m_stat->newGeneration();

	if ((it % JDE_LP)==0) {
	    setStrategyProb();

	    if (m_debug) {
	       printStrategiesStat(m_strategies, strategy_prob);
	    }
	}

	// Select the individual to improve
	i = it % NP;
	ind = m_pop->getInd(i);
	// Get the strategy to use
	int id_strategy = selectStrategy();

        nov = NULL;
	jDE(i, NP, nov, best, m_strategies[id_strategy], it, MaxFES); // calculate i-th offspring

	// function call
	m_new_eval->eval(nov);
	tFitness c = nov->perf();

	// Replace the current individual if it is better
	// offspring has better or equal fitness value
	if (nov->isBetter(ind)) {
	   m_pop->replace(i, nov);

            // Update the best individual
            if (c < best_fit) {
                best_ind = i;
                best_fit = c;
                best = nov->sol();
            }

	    success_memory[id_strategy][it%JDE_LP]+=1;

        }
        else {
	    failure_memory[id_strategy][it%JDE_LP]+=1;
            delete nov;
        }

	if (m_stat)
	  m_stat->endGeneration(best_fit);

	if (NumReduce > 0 && it % (MaxFES / (NumReduce+1)) == (MaxFES/(NumReduce+1) - 1) && NP > 10 && it < MaxFES-1) {
	   // *** new reduction ***
	   m_pop->reduceHalf();
	   NP = m_pop->size();
	}

	// Obtain the better
	best_ind = m_pop->getBest();
	best_fit = m_pop->getInd(best_ind)->perf();

   } // Running

   bestsol = m_pop->getInd(best_ind)->sol();
    
   // Return the best one
   copy(bestsol.begin(), bestsol.end(), sol.begin());
   fitness = best_fit;
   return m_running->numEval();
}

void JDEMC::jDE(const int i, const int NP, tIndividualRealJDEMC* &newsol, tChromosomeReal &best, string strategy, int it, int MaxFES) {
     int r1, r2, r3;
     double tau1 = 0.1;
     double tau2 = 0.1;
     double F, CR;
     int L, n;
     const int  D = m_problem->getDimension();

     tChromosomeReal tmp(m_pop->getInd(i)->sol());

     chooseRs(0,NP,NP,i,r1,r2,r3);     // normal
     tIndividualRealJDEMC *ind = (tIndividualRealJDEMC *) m_pop->getInd(i);
     tIndividualReal *ind1 = m_pop->getInd(r1);
     tIndividualReal *ind2 = m_pop->getInd(r2);
     tIndividualReal *ind3 = m_pop->getInd(r3);
 
     // Each strategy has its own F and CR parameters
     n = mRandomInt(D);
       
     // *** LOWER and UPPER values for strategies *** 
     double F_l = 0.1+sqrt(1.0/NP); // for small NP   
     double CR_l;
     double CR_u = 0.95;

     if (strategy == "jDEbin") { CR_l = 0.00; CR_u=1.0; }
     else if (strategy == "jDEexp" ) { CR_l = 0.3; CR_u=1.0; F_l = 0.5; }
     else { F_l = 0.4; CR_l = 0.7; }  // jDEbest

     if (mRandom()<tau1) {
	F = F_l + mRandom() * (1 - F_l);
     }
     else {
        F = ind->getF(strategy);
     }
     if (mRandom()<tau2) {
        CR = CR_l + mRandom() * (CR_u-CR_l);
     }
     else {
        CR = ind->getCR(strategy);
     }
     assert(strategy == "jDEbin" || strategy == "jDEexp" || strategy == "jDEbest");

     if (strategy == "jDEbin") {
        if (mRandom() < 0.75 && ind2->perf() > ind3->perf())  // COST, sign change 
           F = -F;

        for (L=0; L<D; L++) /* perform D binomial trials */
        {
           if ((mRandom() < CR) || L == (D-1)) /* change at least one parameter */
           {
              tmp[n] = ind1->gen(n) + F*(ind2->gen(n) - ind3->gen(n));
           }
           n = (n+1)%D;
        }
     }
     /*-------DE/rand/1/exp-------------------------------------------------------------------*/
     else if (strategy == "jDEexp") /* strategy DE1 in the techreport */
     {
        if (mRandom() < 0.75 && ind2->perf() > ind3->perf())  // COST, sign change 
           F = -F;

        L = 0;
        //F = 0.5; CR=0.9; 
        do
        {                       
           tmp[n] = ind1->gen(n) + F*(ind2->gen(n)-ind3->gen(n));
           n = (n+1)%D;
           L++;
        }while((mRandom() < CR) && (L < D));
     }
     else if ("jDEbest" == strategy) {         // 
        for (L=0; L<D; L++) /* perform D binomial trials */
        {
           if ((mRandom() < CR) || L == (D-1)) /* change at least one parameter */
           {  
             tmp[n] = (best[n]) + F*(ind2->gen(n)-ind3->gen(n));
           }
           n = (n+1)%D;
        }
     }
    // border check
     m_problem->getDomain()->clip(tmp);
    
     newsol = (tIndividualRealJDEMC*) m_pop->getInstance(tmp);
     newsol->setF(strategy, F);
     newsol->setCR(strategy, CR);
} // end jDE   


int JDEMC::mRandomInt(unsigned size) {
    return m_random->randint(0, size-1);
}

double JDEMC::mRandom() {
    return m_random->rand();
}

void JDEMC::chooseRs(int LB, int UB, int NP, int i, int &r1, int &r2, int& r3) {
   int size = UB - LB; // 
   do { 
       r1 = LB + mRandomInt(size);
   }while(r1==i);
   do { 
       r2 = LB + mRandomInt(size);
       //r2 = mRandomInt(NP);
   }while(r2==i || r2==r1);
   do { 
       r3 = LB + mRandomInt(size);
       //r3 = mRandomInt(NP);
   }while(r3==i || r3==r1 || r3==r2);

}



/**
 * To select the population type
 */
void JDEMC::setPopsize(unsigned int popsize) {
    if (m_pop) {
	delete m_pop;
    }

    m_pop = new PopulationRealJDEMC(m_random,popsize, popsize);
}

int JDEMC::selectStrategy(void)
{
    double r = m_random->rand();
    double aux =0;
    bool end = false;
    int selected=0;

    for(int i=0; i< JDE_NB_OF_STRATEGIES && !end; i++){
        aux += strategy_prob[i];
        if(r <= aux){
             selected = i;
             end  =true;
        }
    }


    return selected;
}

void JDEMC::setStrategyProb(void){
    double sum_of_prob=0;
    double init_value = 1.0/JDE_NB_OF_STRATEGIES;
    int sum_of_failures;
    int sum_of_success;
    
    if (!init_strategies) {
	strategy_prob[0] = init_value;
	strategy_prob[1] = init_value;
	strategy_prob[2] = init_value;
       init_strategies = true;
       return;
    }

    for(int k=0;k<JDE_NB_OF_STRATEGIES;k++)
    {
        sum_of_failures = sum_of_success = 0;
        for(int g=0; g<JDE_LP; g++)
        {
            sum_of_failures += failure_memory[k][g];
            sum_of_success += success_memory[k][g];
        }

	if (sum_of_success == 0 && sum_of_failures == 0) {
	    strategy_prob[k] = 0;
	}
	else {
	    strategy_prob[k]=((double)sum_of_success/(double)(sum_of_failures+sum_of_success))+JDE_EPS;
	}
	sum_of_prob += strategy_prob[k];
    }
    
    for(int k=0;k<JDE_NB_OF_STRATEGIES;k++)
    {
        strategy_prob[k]=strategy_prob[k]/sum_of_prob;
    }


    // Clear the counter
    memset(failure_memory, 0, sizeof(int)*JDE_NB_OF_STRATEGIES*JDE_LP);
    memset(success_memory, 0, sizeof(int)*JDE_NB_OF_STRATEGIES*JDE_LP);
}

