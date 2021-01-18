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

#ifndef _IMUTATION_H 

#define _IMUTATION_H 1

#include "real.h"
#include "domain.h"
#include "random.h"

namespace realea {
/**
 * @class IMutation
 * @ingroup realea_common
 *
 * This class define the mutation type
 *
 * To create a new mutation crossover, the class must inherit from this interface and redefine the 
 * apply method
 */
class IMutation {
	public:
	/**
	 * Make de mutation (changing a gen of the chromosome sol). 
	 *
	 * @param sol the chromosome to update, it is modified 
	 * @param gen the position of gen to change
	 * 
	 */
       virtual tGen mutate(tChromosomeReal &sol, unsigned posi)=0;

	/**
	 * Set the random variable
	 *
	 * It must be speficied if Mutation is going to use it.
	 *
	 * @param random the random generation numbers.
	 */
	void setRandom(Random *random) {
	   m_random = random;
	}
	void setDomain(DomainRealPtr domain) {
	   m_domain = domain;
	}

	virtual ~IMutation(void) {}

protected:
	Random *m_random;
	DomainRealPtr m_domain;
};

}
#endif
