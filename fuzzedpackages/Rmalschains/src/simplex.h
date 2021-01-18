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

#ifndef _SIMPLEX_H

#define _SIMPLEX_H 1

#include "ilocalsearch.h"

namespace realea {

/**
 * @class Simplex
 * apply the LS Simplex 
 * see http://www.grabitech.com/algorithm.htm for details
 */
class Simplex : public ILocalSearch {
    public:
	Simplex(void);
    protected: 

    private:
	unsigned apply(ILSParameters *opt, tChromosomeReal &sol, tFitness &fitness, unsigned itera);
	ILSParameters *getInitOptions(tChromosomeReal &sol);
    private:

    /**
     * Calcalate the centroid
     * Change: m_cum_simplex
     */
    void calculateCentroide(void);

    /**
     * Obtains from m_simplex the extreme values (best and two worst) using the fitness as comparator
     *
     * @param best_index index of best solution
     * @param next_worst_index index of second worst solution
     * @param worst_index index of worst solution
     */
    void getExtremes(ILSParameters *params, int &best_index, int &next_worst_index, int &worst_index);

    /** 
     * Init the p parameter with the first values
     *
     * @param sol initial sol
     * @param fitness sol's fitness
     * @param p LS parameters
     */
    virtual unsigned initParams(tChromosomeReal &sol, tFitness fitness, ILSParameters *param)=0; 

    /**
     * Obtain a new solution using the worst solution and the centroid
     *
     * If the new solution is better than worst solution it is updated.
     *
     * @param worst_point worst solution
     * @param worst_value Fitness of worst solution
     * @param factor Movement factor
     *
     * @return the fitness of new solution
     */
    tFitness move(ILSParameters *p, int posi, double factor); 

    unsigned restart_simplex(ILSParameters *params, int best, unsigned max); 
};

/**
 * @class SimplexDim
 * apply the LS Simplex creating the simplex 
 * see http://www.grabitech.com/algorithm.htm for details
 */
class SimplexDim : public Simplex {
    /** 
     * Init the p parameter with the first values. It creates the
     * simplex structure from the initial point by a mutation in each dimension
     *
     * @param sol initial sol
     * @param fitness sol's fitness
     * @param p LS parameters
     */
    unsigned initParams(tChromosomeReal &sol, tFitness fitness, ILSParameters *param); 
};

/**
 * @class SimplexNeigh
 * apply the LS Simplex creating the simplex from the neighboor
 * see http://www.grabitech.com/algorithm.htm for details
 */
class SimplexNeigh: public Simplex {
    /** 
     * Init the p parameter with the first values. It creates the
     * simplex structure from the initial point by a mutation in each dimension
     *
     * @param sol initial sol
     * @param fitness sol's fitness
     * @param p LS parameters
     */
    unsigned initParams(tChromosomeReal &sol, tFitness fitness, ILSParameters *param); 
};

}
#endif
