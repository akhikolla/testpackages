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

#ifndef _JADE_H

#define _JADE_H 1

#include "ea.h"
#include "distance.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>

# define EPS 0.01
using namespace std;


namespace realea {
/**
 * @class Implementa el método de Diferential Evolution, permitiendo utilizar distintos
 * métodos de cruce (Rand1Exp, Rand2Exp, ...).
 */
class JADE : public ClassEAlgorithm {
   private:
       void cross(PopulationReal *pop, unsigned pos, tChromosomeReal &crom);

       /**
        * Activate debug information (could be very long)
        */
       void setDebug(void);



   public:
       JADE(Random *random);
       unsigned getDefaultPopsize(void) { return 60; }
       ~JADE(void);
       /** 
        * Set the average F (it uses a normal distribution with
        * std=min(0.3, 1-valueF, valueF).
        * 
        * @param meanF
        */
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
       double m_CR;
       double m_F;

       unsigned char m_crossover;
       // Store the meanF
       double m_meanF;
       unsigned m_popReductions;
       double m_meanCR;
       double c;
       double p;
       int m_G;
       unsigned currentPopReduction;
       //deque<tIndividualReal*> m_archive;
       vector<tIndividualReal*> m_archive;
};

}

#endif
