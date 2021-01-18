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

#ifndef _CMAESHANSEN_H 

#define _CMAESHANSEN_H 1

#include "ilocalsearch.h"

namespace realea {

/**
 * @class CMAESHansen
 *
 * @brief Local Search method CMAESHansen defined by Hansen. 
 *
 * It has been reimplementated by Daniel Molina, from the implementation in  Matlab from Hansen
 */
class CMAESHansen : public ILocalSearch {
    public:
	/**
	 * Constructor
	 *
	 * @param fconfig file config
	 */
	CMAESHansen(string fconfig);
	/**
	 * Search in the total range
	 *
	 * @param factor ratio of search
	 */
	void searchRange(double factor);
	/**
	 * Enable the local search bound checking
         * 
	 */
	void enableBoundChecking(void);

        /**
         * Select my current random seed when the LS is applied
         */
        void setMyRandom(void);
	/**
	 * Disable the local search bound checking
	 */
	void disableBoundChecking(void);
	/**
	 * Search in the population neighbood
	 *
	 * @param factor ratio of search
	 */
	void searchNeighborhood(double ratio);
        /**
         * Set the lamba, to use (0 to default)
         */
        void setPopsize(int lambda);
        /**
         * Set the number of parents (mu)
         */
        void setParentsSize(int mu);

        int getPopsize();
        int getParentsSize();

	unsigned apply(ILSParameters *opt, tChromosomeReal &sol, tFitness &fitness, unsigned itera);
	ILSParameters *getInitOptions(tChromosomeReal &sol);

    private:
        int m_lambda;
        int m_mu;
	double m_rfactor;
	double m_nfactor;
	string m_fconfig;
	bool m_debug;
        bool m_setmyrandom;
	string m_bound_strategy;
};

}
#endif
