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

#ifndef _MTSI_LS_H 

#define _MTSI_LS_H 1

#include "ilocalsearch.h"

namespace realea {

/**
 * This LS method is a variant of the LS from MTS algorithm with a fixed gen explorations
 *
 * This change dimension to dimension
 */
class MTSILS : public IParallelLocalSearch {
    public:
	MTSILS(double maxdelta, double mindelta, unsigned maxgen);
    private:
	unsigned apply(ILSParameters *opt, tChromosomeReal &sol, tFitness &fitness, unsigned itera);
	ILSParameters *getInitOptions(tChromosomeReal &sol);
	ILSParameters *recoverOptions(tGen *params, unsigned size);

    /**
     * Store in a real vector the values of a parameters. 
     *
     * @param params parameter to store
     * @param paparams reference to sequencial parameters (can be NULL)
     * @param size reference to size (always is right, even if params == NULL)
     *
     * @see recoverOptions
     */
    void storeOptions(ILSParameters *params, tGen **paparams, unsigned *psize);

    private:
	double m_maxdelta;
	double m_mindelta;
	unsigned m_maxgen;
};

}
#endif
