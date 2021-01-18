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

#ifndef _SRANDOM_H

#define _SRANDOM_H

#include "random.h"

/**
 * @class ORandom
 *
 * @ingroup realea_common
 *
 * @brief This class generates number randomly using a pseudogenerator that uses a PRIME multiplication
 */

class ORandom : public IRealRandom {
       /**
	* init the seed 
	*
	* @param seed (value != 0)
	*/
       void setSeed(unsigned long seed);

       /**
	* @return A random double between 0 and 1 
	*/
       virtual double rand(void);

      /**
	* @return the actual seed
	*/
       unsigned long getSeed(void);
private:
	unsigned long m_seed; /*< seed */
};

#endif
