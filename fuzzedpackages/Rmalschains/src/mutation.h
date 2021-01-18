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

#ifndef _MUTATION_H 

#define _MUTATION_H 1

#include "imutation.h"

namespace realea {
/**
 * @class Mutation, mutación operator (for AGs). 
 *
 * @ingroup realea_common
 *
 * Usage:
 * 
 * m_mutation->apply(sol);
 *
 */
class Mutation {
   public:
       /**
	* Constructor.
	*
	* @param mut IMutation to be applied if it is decided
	*
	* @param ratio Ratio of individuals to be applied
	*/
       Mutation(IMutation *mut, double ratio=0.125);

	/**
	 * Set the random variable
	 *
	 * It must be speficied if Mutation is going to use it.
	 *
	 * @param random the random generation numbers.
	 */

       void setRandom(Random *random);
       void setDomain(DomainRealPtr domain);
       /**
        * Set a ratio to the mutation
        *
        * @param ratio new ratio (between 0 and 1)
        */
       void setRatio(double ratio);

       virtual ~Mutation(void);
       /**
	* This method check if must be applied the chromosome, in that case it is updated.
	*
	* @return true if the solution has been changed.
	*/
       bool apply(tChromosomeReal &sol);
    private:
       double m_ratio;
       IMutation *m_mutation;
       DomainRealPtr m_domain;

    protected:
       Random *m_random;
};

typedef Mutation* MutationPtr;

class MutationBGA : public IMutation {
    public:
       virtual tGen mutate(tChromosomeReal &sol, unsigned pos);
};

}
#endif
