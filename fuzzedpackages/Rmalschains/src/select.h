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

#ifndef _SELECT_H

#define _SELECT_H 1

#include "populationreal.h"
#include "random.h"
#include "signal.h"

namespace realea {
/**
 * Esta clase permite seleccionar los padres a cruzar (para el AGE)
 *
 * @class ISelect
 * @ingroup realea_ea_ssga
 *
 * @brief abstract the parent selection method
 */
class ISelect : public IReset {
public:
    /**
     * Select two new parents from a population
     *
     * @param pop population from select the individuals
     *
     * @param mom mother of the crossover
     *
     * @param dad father of the crossover
     *
     */
    virtual void select(PopulationRealPtr pop, unsigned *mom, unsigned *dad)=0;
    /**
     * Set the random number generators
     *
     * @param random random generator
     */
     void setRandom(Random *random) {
	  m_random = random;
     }
     void setDomain(DomainRealPtr domain) {
	 m_domain = domain;
     }
protected:
    Random *m_random;
    DomainRealPtr m_domain;
};

/**
 * @class SelectNAM
 *
 * @brief Negative Assortative Mating parent mechanism
 *
 * It select the first one (mother) randomly, for the election of the father, it select randomly
 * nam other individuals and choses the more diferent to the first one (using the Euclidean distance)
 */
class SelectNAM : public ISelect {
public:
    /**
     * Constructor.
     *
     * @param random Random number generator
     *
     * @param nam number of individuals considered by the selection of the second individual
     */
    SelectNAM(unsigned nam);
    virtual void select(PopulationRealPtr pop, unsigned *mom, unsigned *dad);
private:
    unsigned m_num;
};

/**
 * @class SelectTournament
 *
 * @brief Tournament parent mechanism
 *
 * It select the two parents using a 2 binary tournament
 */
class SelectTournament : public ISelect {
public:
    /**
     * Constructor.
     *
     * @param random Random number generator
     *
     * @param nam number of individuals in each tournament
     */
    SelectTournament(unsigned nam);
    virtual void select(PopulationRealPtr pop, unsigned *mom, unsigned *dad);
private:
    unsigned m_num;
};

}
#endif
