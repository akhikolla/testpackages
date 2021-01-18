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

#ifndef _SADEAF_H

#define _SADEAF_H 1

#include "ea.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>

#ifndef NB_OF_STRATEGIES 
#define NB_OF_STRATEGIES 4
#define EPS 0.01
#define LP 50
#endif

using namespace std;


namespace realea {
/**
 * @class Implementa el método de Diferential Evolution, permitiendo utilizar distintos
 * métodos de cruce (Rand1Exp, Rand2Exp, ...).
 */
class SADEAF : public ClassEAlgorithm {
   private:
       void cross(PopulationReal *pop, unsigned pos, tChromosomeReal &crom, int strategy);
       void crossRand1Bin(PopulationReal *pop, unsigned pos, tChromosomeReal &crom);
       void crossRandToBest2Bin(PopulationReal *pop, unsigned pos, tChromosomeReal &crom);
       void crossRand2Bin(PopulationReal *pop, unsigned pos, tChromosomeReal &crom);
       void crossCurrentToRand1(PopulationReal *pop, unsigned pos, tChromosomeReal &crom);

       void setStrategyProb();
       int selectStrategy();
       /**
        * Activate debug information (could be very long)
        */
       void setDebug(void);

       void printStrategyProb();
       void printCRmk();
       void printCRk();
       void printSuccessMemory();
       void printFailureMemory();
       void printCRMemory();



   public:
       SADEAF(Random *random);
       unsigned getDefaultPopsize(void) { return 60; }
       ~SADEAF(void);
       /** 
        * Set the average F (it uses a normal distribution with
        * std=min(0.3, 1-valueF, valueF).
        * 
        * @param meanF
        */    
       void setAverageF(double meanF);
       unsigned init(void);
       unsigned realApply(tChromosomeReal &sol, tFitness &fitness);
       void setCross(ICrossBinaryPtr cross);

       /** 
        * Set the population reductions (Default it is disable). 
        * 
        * @param num reductions number during the evaluations
        */
	void setPopReductions(unsigned num);
       //void setLP(int l);
    private:

       //int LP;
       unsigned char m_crossover;
       // Store the meanF
       double initial_meanF;
       double m_meanF;
       unsigned m_popReductions;
       double m_CR;
       double m_F;
       double m_K;
       int m_G;
       unsigned currentPopReduction;

       int failure_memory[NB_OF_STRATEGIES][LP];

       int success_memory[NB_OF_STRATEGIES][LP];

       double strategy_prob[NB_OF_STRATEGIES];
       
       double CR_memory[NB_OF_STRATEGIES][LP];
       double CRk[NB_OF_STRATEGIES];
       double CRmk[NB_OF_STRATEGIES];

       double F_memory[LP];
};

}

#endif
