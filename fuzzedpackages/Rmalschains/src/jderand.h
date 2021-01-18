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

#ifndef _JDERAND_H

#define _JDERAND_H 1

#include "ea.h"
#include "populationjderand.h"

namespace realea {
/**
 * @class Implements the jDEdynNP-FCR, proposed by Janez Brest et al in:
 *
 * Janez Brest; Ales Zamuda; Iztok Fister; Mirjam Sepesy Maucec. Large Scale Global Optimization using Self-adaptive 
 * Differential Evolution Algorithm. Evolutionary Computation, 2010. WCCI 2010 IEEE World Congress on Computational Intelligence. IEEE Congress on
July, 18-23, 2010 Page(s): 3097 - 3104
 */
class JDERand : public ClassEAlgorithm {
   private:
       /**
	* Obtain the mutation strategy to use in function of the current iteration
	*
	* @param iteration current iteration
	* @param maxFES maximum iteration
	* @return string with the mutation to apply to the individual
	*/
   public:
       JDERand(Random *random, int reduced=0);
       unsigned getDefaultPopsize(void) { return 100; }
       ~JDERand(void);
       unsigned init(void);
       unsigned realApply(tChromosomeReal &sol, tFitness &fitness);
       void setPopsize(unsigned int popsize);

       /**
	* Allow to set the strategy
	*/
       void setStrategy(string name);
       
    private:
        /**
         * Obtain the current strategy to apply
         *
         * @param NP popsize
         * @param iteration current iteration
         * @param MaxFEs maximum FEs
         * 
         * @return the strategy name to apply in the current state of the run
         */
        string getStrategy();

        /**
         *
         * @param i individual position
         * @param NP popsize
         * @param newsol new solution
         * @param best current best solution
         * @param strategy to apply
         * @param it iteration number
         * @param MaxFES maximum FE
         */
        void jDE(const int i, const int NP, tIndividualRealJDERand* &newsol, tChromosomeReal &best, string strategy, int it, int MaxFES);
private:
        /**
         * Get a random value between [0, size-1]
         * @param size
         * @return random value
         */
        int mRandomInt(unsigned size);
        /**
         *
         * @return a random value between [0, 1]
         */
        double mRandom(void);
        /**
         * Choose the random values
         * @param LB initial value
         * @param UB final value
         * @param NP popsize
         * @param i current individual pos
         * @param r1 new individidual
         * @param r2 new individidual
         * @param r3 new individidual
         */

        void chooseRs(int LB, int UB, int NP, int i, int &r1, int &r2, int& r3);

       int NumReduce; // 0.. no reduction, 1, ... max. num. of reductions
       int extDim; // extended dim to store control parameters
       string m_strategy;

};

}

#endif
