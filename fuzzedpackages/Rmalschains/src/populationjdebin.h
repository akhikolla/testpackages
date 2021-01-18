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

#ifndef _INDIVIDUALJDEBin_H

#define _INDIVIDUALJDEBin_H 1

#include "individual.h"
#include "domain.h"
#include "populationreal.h"
#include <bitset>

#define BITS_GEN 30

namespace realea {

class tIndividualRealJDEBin : public tIndividualReal {
public:
    tIndividualRealJDEBin (const tChromosomeReal &com);
    tIndividualRealJDEBin (const tChromosomeReal &com, double fitness);
    virtual ~tIndividualRealJDEBin(void);
    double getF(string strategy);
    double getCR(string strategy);
    void setF(string strategy, double F);
    void setCR(string strategy, double CR);

private:
    map<string,double> m_F; // self-ad. F init. for each strategy
    map<string,double> m_CR; // self-ad CR init. for each strategy
};


/**
 * This class allow us to define a population with the
 * binary Codification
 */
class PopulationRealJDEBin : public PopulationReal {
   public:
    PopulationRealJDEBin(Random *random,unsigned int max, unsigned int pob);
    void restart(DomainRealPtr domain);

   private:
       tIndividualReal* getInstance(tChromosomeReal &crom);
       tIndividualReal* getInstance(tChromosomeReal &crom, double fitness);
};

}
#endif
