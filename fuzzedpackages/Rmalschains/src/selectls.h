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

#ifndef _SELECT_IMPROVEMENT_LS_H 

#define _SELECT_IMPROVEMENT_LS_H 1

#include "populationreal.h"

namespace realea {

class SelectImprovementLS {
public:
	/**
	 * Select the individuals of the population that could be improved by the LS method
	 *
	 * @param pop population to search
	 *
	 * @param subpop an deque of individuals with the elements to improve
	 */
    virtual void getIndsToImprove(PopulationReal *pop, deque<tIndividualReal*> &subpop)=0;

	/**
	 * Return the individual id to improve 
	 *
	 * @param individuals individuals to consider
	 *
	 * @return identification of the chosen individual 
	 *
	 */
    virtual unsigned selectIndToImprove(deque<tIndividualReal*> &individuals)=0; 
    virtual ~SelectImprovementLS(void) {}
};

}

#endif
