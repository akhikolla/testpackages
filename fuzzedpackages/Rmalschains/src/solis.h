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

#ifndef _SOLIS_H 

#define _SOLIS_H 1

#include "ilocalsearch.h"

namespace realea {


class SolisWets : public IParallelLocalSearch {
    public:SolisWets(void);
    void setDelta(double maxdelta);
    void setDelta(double mindelta, double maxdelta);
    protected: 
    /**
     * @var Máximo número de mejoras
     */
    static const unsigned maxSuccess=5;

    /**
     * @var Mínimo número de mejoras consecutivas
     */
    static const unsigned maxFailed=3;

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
	/**
	 * Obtain a neighbour solution from a solution given.
	 *
	 * @param p Solis Wets parameters
	 * @param sol current solution
	 * @param dif  calculate ddifferences vector, output
	 * @param newsol new solution, output
	 *
	 * @return fitness of newsol
	 */
	tFitness getNeighbour(ILSParameters *param, tChromosomeReal &sol, tChromosomeReal &dif, tChromosomeReal &newsol);

	double m_maxdelta;
	double m_mindelta;
};

}
#endif
