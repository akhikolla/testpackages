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

#ifndef _DE_H

#define _DE_H 1

#include "ea.h"

namespace realea {
/**
 * @class Implementa el método de Diferential Evolution, permitiendo utilizar distintos
 * métodos de cruce (Rand1Exp, Rand2Exp, ...).
 */
class DE : public ClassEAlgorithm {
   private:
       void crossBin(PopulationReal *pop, unsigned pos, tChromosomeReal &crom);
       void crossExp(PopulationReal *pop, unsigned pos, tChromosomeReal &crom);
       void cross(PopulationReal *pop, unsigned pos, tChromosomeReal &crom);

   public:
       DE(Random *random);
       unsigned getDefaultPopsize(void) { return 60; }
       ~DE(void);
       unsigned init(void);
       unsigned realApply(tChromosomeReal &sol, tFitness &fitness);
       void setCross(ICrossBinaryPtr cross);
       void setCrossoverBin(void);
       void setCrossoverExp(void);
       void setF(double F);
       void setCR(double CR);
    private:
       unsigned char m_crossover;
       double m_CR;
       double m_F;
};

}

#endif
