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

#ifndef _CMAES_H 

#define _CMAES_H 1

#include "ilocalsearch.h"
#include "newutil.h"

namespace realea {

/**
 * @class CMAES
 *
 * @brief Local Search method CMAES defined by Hansen. 
 *
 * It has been reimplementated by Daniel Molina, from the implementation in  Matlab from Hansen
 */
class CMAES : public ILocalSearch {
    public:
	/**
	 * Constructor
	 *
	 * @param factor (reduction factor from min distance in population)
	 */
	CMAES(void);
	/**
	 * Search in the total range
	 *
	 * @param factor ratio of search
	 */
	void searchRange(double factor);
	/**
	 * Search in the population neighbood
	 *
	 * @param factor ratio of search
	 */
	void searchNeighborhood(double ratio);
	unsigned apply(ILSParameters *opt, tChromosomeReal &sol, double &fitness, unsigned itera);
	ILSParameters *getInitOptions(tChromosomeReal &sol);

    private:
	double m_rfactor;
	double m_nfactor;
	bool m_debug;
};

}
#endif
