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

#ifndef _UTIL_H

#define _UTIL_H 1

#include "define.h"
#include "populationreal.h"

using realea::tChromosomeReal; 
using realea::PopulationReal; 


//double distance(tGen *x, tGen *y, unsigned size);
double distreal(const realea::tChromosomeReal &x, const realea::tChromosomeReal &y, bool *checkGen=NULL);

double distanceMin(const realea::tChromosomeReal &x, realea::PopulationReal *pop, unsigned *posmin);
double distanceMax(const realea::tChromosomeReal &x, realea::PopulationReal *pop, unsigned *posmin);

unsigned getNeigh(const tChromosomeReal &x, PopulationReal *pop, double min); 

/**
 * Obtain the vector distance from an individual to the rest
 *
 * @param x individual
 * @param pop population
 * @param mindist output with the mean difference distance
 */
void min_vector_distance(const tChromosomeReal &x, PopulationReal *pop, vector<tReal> &mindist); 

/**
 * This method return the individual closer in each solution (ignoring the
 * equal dimension values
 *
 * @param x initial solution
 * @param pop population
 * @param mindist position distance for each dimension
 */
void min_dim_distance(const tChromosomeReal &sol, PopulationReal *pop, vector<unsigned> &minind);

#endif
