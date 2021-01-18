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

#ifndef _DOMAIN_H 

#define _DOMAIN_H 1

#include "real.h"
#include "random.h"

namespace realea {
/**
 * @class DomainReal
 * @ingroup realea_common
 *
 * Represents the search domain. 
 *
 * It stores the maximum and minimum values for each dimension.
 */
class DomainReal {
   public:
       /**
	* Constructor without specify the max and min values.
	* If this constructor is applied, setValues must be specify for each value
	*
	* @param dim dimension of chromosomes
	* @see setValues
	*/
	DomainReal(unsigned dim);
	~DomainReal(void);

	/**
	 * Specifies the min and max value for a gen
	 *
	 * @param gen gen to configurate (from 0 to ndim() -1)
	 * @param min minimum value for this gen
	 * @param max maximum value for this gen
	 * @param check if it must be checked the position gen
	 * @see getValues
	 */
	void setValues(unsigned int gen, tReal min, tReal max,bool check=true);

	/**
	 * Set the dimensions to search
	 *
	 * @param searchDim vector of boolean indicating the gens to be searched
	 *
	 * @param dim dimension of the vector
	 */
	void setSearchDomain(bool *searchDim, unsigned dim);

	/**
	 * Get the dimension to search
	 *
	 */
	void getSearchDomain(bool *searchDim, unsigned dim);

	/**
	 *  Recovers the min and max value for a gen
	 *  
	 *  @param gen gen to consult (from 0 to ndim() -1)
	 *  @param pmin minimun value of this gen, output
	 *  @param pmax maximun value of this gen, output
	 *  @param check if it must be checked the position gen
	 */
	void getValues(unsigned int gen, tReal *pmin, tReal *pmax, bool check=true);

	/**
	 * Return if the dimension index can be changed during the search
	 *
	 * @param dim dimension index (0..ndim)
	 *
	 * @return true if the current dimension dim must be used into the search (false in othercase)
	 */
	bool canBeChanged(unsigned dim);

	/**
	 * clip the real (checking that it is between the [min,max] values for this gen)
	 *
	 * @param gen gen to check
	 * @param value value to check
	 *
	 * @return the min if (value < min), value if (value <= max); max if (value > max), original
	 * value othercase
	 *
	 */
	tReal clip(unsigned int gen, tReal value, bool check=true);

	/**
	 * Set the scale of the domain search around a solution
	 *
	 * @param center center of the new domain. 
	 * @param scale scale of the domain search. 1->not changed, 0.1->10% of previous scale. 
	 *
	 * Never it allow to over the original bounds
	 *
	 */
	void setDomainCenter(tChromosomeReal center, double scale);

	/**
	 * Clip each one of the gens of the indicated chromosome
	 *
	 * @param crom chromosome to check, it can be modified
	 */
	void clip(tChromosomeReal &crom);
        void clip(double *crom);

	/**
	 * check if the chromosome folow the domain restrictions
	 *
	 * @param crom chrosomome to check, it can't be changed
	 * @return true if the chromosome has valid values
	 */
	bool check(const tChromosomeReal &crom);
	
	/**
	 * Can obtain a new chromosome
	 *
	 * @var random generator of random values
	 * @var crom chromosome to generate
	 */
	void getInit(Random *random, tChromosomeReal &crom);
	
	/**
	 * @return the dimensionality
	 */
	unsigned getDimension(void) {
	    return m_dim;
	}

	/**
	 * Set the bound. It is set the clip check the solutions
	 * 
	 */
	void setBounds(void) {
	     m_isbound = true;		
	}
	
	void setNotBounds(void) {
	     m_isbound = false;	
	}

	bool isBound() {
	     return m_isbound;
	}

   private:
	/**
	 * Can obtain a chromosome randomly
	 *
	 * @var random generator of random values
	 * @var crom chromosome to be generated
	 */
	void getInitRandom(Random *random, tChromosomeReal &crom); 
	
	/**
	 * check if the gen is value, using asserts
	 */
	void checkGen(unsigned int gen);
	tChromosomeReal m_mins; /// min values vector
	tChromosomeReal m_maxs; /// max values vector
	unsigned int m_dim; /// dimensionality
	bool m_isbound; /// if true must be checked the bounds
	unsigned m_search_ini; // The initial dimension search 
	unsigned m_search_fin; // The final dimension search

	bool *m_check_dim; // The dimension search
};

typedef DomainReal* DomainRealPtr;

}
#endif
