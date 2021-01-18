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

#include "distance.h"
#include <cmath>
#include <cassert>

using namespace realea;

//double distance(tGen *x, tGen *y, unsigned size) {
//    double dist = 0;
//
//    for (i = 0; i < size; i++) {
//	dist += (x[i]*y[i])*(x[i]*y[i]);
//    }
//
//    return sqrt(dist);
//}

double distreal(const tChromosomeReal &x, const tChromosomeReal &y, bool *checkGen) {
    double dist = 0;
    unsigned size = x.size();
    assert(x.size() == y.size());

    for (unsigned i = 0; i < size; i++) {
        if (checkGen == NULL || checkGen[i]) {
	    dist += (x[i]-y[i])*(x[i]-y[i]);
	}
    }

    return sqrt(dist)/size;
}

double distanceMin(const tChromosomeReal &x, PopulationReal *pop, unsigned *posmin) {
   tIndividualRealPtr ind;
   double dist, lowest;
   unsigned i;

    if (pop->size() == 0) {
       throw new string("dist:Error, popsize is zero");
    }

    *posmin = 0;
    lowest = 0;

   for (i = 0; i < pop->size(); i++) {
	ind = pop->getInd(i);
	dist = distreal(x, ind->sol());

	if (dist > 0 && (lowest == 0 || dist < lowest) ) {
	   *posmin = i;
	   lowest = dist;     
	}
   }

   return lowest;
}

double distanceMax(const tChromosomeReal &x, PopulationReal *pop, unsigned *posmin) {
   tIndividualRealPtr ind;
   double dist, greatest;
   unsigned i;

    if (pop->size() == 0) {
       throw new string("dist:Error, popsize is zero");
    }

    ind = pop->getInd(0);
    greatest = distreal(x, ind->sol());

   for (i = 0; i < pop->size(); i++) {
	ind = pop->getInd(i);
	dist = distreal(x, ind->sol());

	if (dist > greatest) {
	   *posmin = i;
	   greatest = dist;     
	}
   }

   return greatest;
}


void vector_distance(const tChromosomeReal &x, const tChromosomeReal &y, vector<tGen> &mindist) {
    unsigned size = x.size();
    assert(x.size() == y.size());

    for (unsigned i = 0; i < size; i++) {
	mindist[i] = fabs(x[i]-y[i]);
    }
}


void min_vector_distance(const tChromosomeReal &x, const tChromosomeReal &y, vector<tGen> &mindist) {
    double dist = 0;
    unsigned size = x.size();
    assert(x.size() == y.size());

    for (unsigned i = 0; i < size; i++) {
	dist = fabs(x[i]-y[i]);
	
	if (dist < mindist[i] && dist > 0) {
	   mindist[i] = dist;
	}
    }

}


unsigned getNeigh(const tChromosomeReal &x, PopulationReal *pop, double min) {
   double dist, distMin;
   unsigned posmin,i,size;
   posmin = 0;

    distMin = -1;
    size = pop->size();

   for (i = 0; i < size; i++) {
	dist = distreal(x, pop->getInd(i)->sol()); 

	if (dist > min && (distMin < 0 || dist < distMin)) {
	   posmin = i;
	   distMin = dist;
	}
   }

   assert(distMin > 0);
   return posmin;
}

void min_vector_distance(const tChromosomeReal &x, PopulationReal *pop, vector<tGen> &mindist) {
   tIndividualRealPtr ind;
   vector<double> distMin(x.size());
   unsigned i, initial;

    if (pop->size() == 0) {
       throw new string("dist:Error, popsize is zero");
    }

    initial = 0;
    ind = pop->getInd(initial);

    if (ind->sol() != x) {
	initial += 1;
	ind = pop->getInd(initial);
    }

    vector_distance(x, ind->sol(), mindist);

   for (i = initial+1; i < pop->size(); i++) {
	ind = pop->getInd(i);
	min_vector_distance(x, ind->sol(), mindist);
   }

}

void min_dim_distance(const tChromosomeReal &sol, PopulationReal *pop, vector<unsigned> &minind) {
    tChromosomeReal sol_ind;
    unsigned i, dim;
    unsigned ndim=sol.size();
    double dist;

    vector<double> mindist(ndim);
    assert(minind.size()==ndim);
    fill(mindist.begin(), mindist.end(), 0);

    if (pop->size() == 0) {
       throw new string("dist:Error, popsize is zero");
    }

    for (i = 0; i < pop->size(); i++) {
	sol_ind = pop->getInd(i)->sol();
	
	for (dim = 0; dim < ndim; ++dim) {
	    dist = fabs(sol_ind[dim]-sol[dim]);
	    
	    if (mindist[dim] == 0 || (dist > 0 && (dist < mindist[dim]))) {
		mindist[dim] = dist;
		minind[dim] = i;
	    }
	}
	
   }
    
}
