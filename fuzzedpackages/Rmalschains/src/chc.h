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

#ifndef _CHC_H 

#define _CHC_H 1

#include "populationchc.h"
#include "cross.h"
#include "iea.h"

namespace realea {
/**
 * @class Implementa el método CHC
 *
 */
class CHC : public ICrossEAlgorithm {
public:
    CHC(Random *random);
    ~CHC(void);
    unsigned getDefaultPopsize(void);

    unsigned init(void);
    unsigned realApply(tChromosomeReal &sol, tFitness &fitness);
    void setPopsize(unsigned popsize);

private:
    unsigned cross(unsigned threshold);
    unsigned int distHamming(tIndividualRealPtr ind1, tIndividualRealPtr ind2);
    int m_init_threshold, m_threshold; 
};

class CHCShow : public Statistics {
public:
    CHCShow(void) : Statistics(1) {
       activeEvent("Restart");
    }
    CHCShow(unsigned num) : Statistics(num) {
       activeEvent("Restart");
    }

};

}

#endif
