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

#include "jderand.h"
#include "random.h"
#include "populationreal.h"
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <cstdio>

#define CopyVector(a,b) memcpy((a),(b),nDim*sizeof(double))
#define Element(a,b,c)  (a->getInd(b).gen(c))

using namespace std;
using namespace realea;

JDERand::JDERand(Random *random, int reduced) : ClassEAlgorithm(random) {
    NumReduce = reduced;    // 7.2.2010 pmax := NumReduce+1 ==> pmax = 3+1 = 4
    m_strategy = "jDEbin";
}

JDERand::~JDERand(void) {
}


unsigned JDERand::init(void) {
    // init/check the parameters
    int dim = m_problem->getDimension();
    extDim = dim + 30; // extended dim to store control paramaters

    // Init the population
    m_pop->reset(m_problem->getDomain());
    // Init the different elements: Running, Cross ...
    reset();
    // Eval the population
    m_pop->eval(m_init_eval);
    return m_running->numEval();
}

unsigned JDERand::realApply(tChromosomeReal &sol, tFitness &fitness) {
    const int D=m_problem->getDimension();
    tChromosomeReal crom(D);
    tChromosomeReal best(D);
    tIndividualRealJDERand* nov;
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

	// Select the individual to improve
	i = it % NP;
	ind = m_pop->getInd(i);


	// Get the strategy to use
	string strategy = getStrategy();

        nov = NULL;
	jDE(i, NP, nov, best, strategy, it, MaxFES); // calculate i-th offspring
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

        }
        else {
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

string JDERand::getStrategy() {
    string strategy;

    double  draw = m_random->rand();

    if (draw < (1/3))
      strategy = "jDEbest";
    else if (draw < (2/3))
      strategy = "jDEbin";
    else
        strategy = "jDEexp";

    return strategy;
}

void JDERand::jDE(const int i, const int NP, tIndividualRealJDERand* &newsol, tChromosomeReal &best, string strategy, int it, int MaxFES) {
     int r1, r2, r3;
     double tau1 = 0.1;
     double tau2 = 0.1;
     double F, CR;
     int L, n;
     const int  D = m_problem->getDimension();

     tChromosomeReal tmp(m_pop->getInd(i)->sol());

     chooseRs(0,NP,NP,i,r1,r2,r3);     // normal
     tIndividualRealJDERand *ind = (tIndividualRealJDERand *) m_pop->getInd(i);
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
     /*-------DE/current-to.best/1/bin-------------------------------------------------------------------*/
     else if ("jDEbest" == strategy) {         // 
        for (L=0; L<D; L++) /* perform D binomial trials */
        {
           if ((mRandom() < CR) || L == (D-1)) /* change at least one parameter */
           {  
             tmp[n] = ind->gen(n)+F*(best[n]-ind->gen(n)) + F*(ind2->gen(n)-ind3->gen(n));
           }
           n = (n+1)%D;
        }
     }
    // border check
     m_problem->getDomain()->clip(tmp);
    
     newsol = (tIndividualRealJDERand*) m_pop->getInstance(tmp);
     newsol->setF(strategy, F);
     newsol->setCR(strategy, CR);
} // end jDE   


int JDERand::mRandomInt(unsigned size) {
    return m_random->randint(0, size-1);
}

double JDERand::mRandom() {
    return m_random->rand();
}

void JDERand::chooseRs(int LB, int UB, int NP, int i, int &r1, int &r2, int& r3) {
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
void JDERand::setPopsize(unsigned int popsize) {
    if (m_pop) {
	delete m_pop;
    }

    m_pop = new PopulationRealJDERand(m_random,popsize, popsize);
}
