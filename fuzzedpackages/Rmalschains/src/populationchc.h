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

#ifndef _INDIVIDUALCHC_H

#define _INDIVIDUALCHC_H 1

#include "individual.h"
#include "domain.h"
#include "populationreal.h"
#include <bitset>

#define BITS_GEN 30

namespace realea {

class tIndividualRealCHC : public tIndividualReal {
public:
    tIndividualRealCHC (const tChromosomeReal &com);
    tIndividualRealCHC (const tChromosomeReal &com, double fitness);
    virtual ~tIndividualRealCHC(void);
    /**
     * This method calculates the gray representation
     */
    void calculateBin(DomainRealPtr domain);


   /**
     * Calcula la distancia de Hamming entre este elemento y
     * el indicado
     */
    unsigned distHamming(tIndividualRealCHC *other);
    unsigned distHammingOpt(tIndividualRealCHC *other);
private:

    char *getBin(void);
    unsigned getBinSize(void);
 
    char *m_codbin;
    vector< bitset<BITS_GEN> > m_codbin_opt;
    unsigned int m_codbin_size;
    bool m_is_codbin;
};


/**
 * This class allow us to define a population with the
 * binary Codification
 */
class PopulationRealCHC : public PopulationReal {
   public:
    PopulationRealCHC(Random *random,unsigned int max, unsigned int pob);
    void restart(DomainRealPtr domain);

   private:
       tIndividualReal* getInstance(tChromosomeReal &crom);
       tIndividualReal* getInstance(tChromosomeReal &crom, double fitness);
};

}
#endif
