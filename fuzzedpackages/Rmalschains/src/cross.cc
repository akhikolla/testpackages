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

#include "cross.h"
#include "problem.h"
#include "random.h"
#include <cassert>
#include <cmath>

using namespace realea;
using namespace realea::internal;

void CrossBinary::operator()(tIndividualReal *mom, tIndividualReal *dad, tChromosomeReal &child) {
    mom->incremCount("cross");
    dad->incremCount("cross");
    operator()(mom->sol(), mom->perf(), dad->sol(), dad->perf(), child);
}

void CrossBinary::operator()(const tChromosomeReal &mom, tFitness fit_mom, const tChromosomeReal &dad, tFitness fit_dad, tChromosomeReal &children) {
    (*m_cross)(mom,fit_mom, dad, fit_dad, children);
}

CrossBinary::CrossBinary(ICrossBinaryPtr cross) : m_cross(cross) {
}

CrossBinary::~CrossBinary(void) {
	delete m_cross;
}
 
 
CrossBLX::CrossBLX(double alpha) : m_alpha(alpha) {
   assert(alpha > 0);
}

CrossBLX::~CrossBLX(void) {
}



void CrossBLX::operator()(const tChromosomeReal &mom, tFitness fit_mom, const tChromosomeReal &dad, tFitness fit_dad, tChromosomeReal &children) {
    double puntox, puntoy, I, A1, B1;
    double rate;
    unsigned i, size = dad.size();
    tGen min, max;

   assert(dad.size() == mom.size());
// First copy all gen from mother
   copy(mom.begin(), mom.end(), children.begin());

   if (fit_mom == fit_dad) {
      rate = 0.5;
   }
   else {
      rate = fit_dad/(fit_mom+fit_dad);
   }

   rate = 0.5;

    /* REALIZAR EL CRUCE PARA Pm*Tp */
   for (i=0; i< size; i++)
     {
       if (!m_domain->canBeChanged(i)) {
	  if (m_random->rand() <= rate) {
	     children[i] = mom[i];
	  }
	  else {
	     children[i] = dad[i];
	  }

       }

       m_domain->getValues(i, &min, &max);
       puntox=mom[i];  puntoy=dad[i];
       if (puntox > puntoy) {
	  swap(puntox,puntoy); 
       }

       I=(puntoy-puntox)*m_alpha;
       assert(I >= 0);

       A1 = puntox-I;

       if (A1 < min && m_domain->isBound())
	  A1 = min;

       B1 = puntoy+I;

	if (B1 > max && m_domain->isBound()) {
	    B1 = max;
	}

       children[i]=A1 + m_random->rand()*(B1-A1);
     }

}

CrossDim::CrossDim(double alpha, double p_r) : m_alpha(alpha), m_pr(p_r) {
   assert(alpha > 0);
}

CrossDim::~CrossDim(void) {
}



void CrossDim::operator()(const tChromosomeReal &mom, tFitness fit_mom, const tChromosomeReal &dad, tFitness fit_dad, tChromosomeReal &children) {
    double puntox, puntoy, I, A1, B1;
    double rate;
    unsigned i, size = dad.size();
    vector<bool> changed(size);
    tGen min, max;

   assert(dad.size() == mom.size());
// First copy all gen from mother
   copy(mom.begin(), mom.end(), children.begin());

   if (m_random->rand() < 0.5) {
      rate = 0;
   }
   else {
      rate = 1;
   }

   fill(changed.begin(), changed.end(), false);
   unsigned int n = m_random->randint(0, size-1);

//   for (i = 0; i < n && m_random->rand() <= m_pr; i++) {

     for (i = 0; i < n; i++) {
	if (m_random->rand() < m_pr) {
	   changed[(i+n)%size] = true;
	}
	
     }

    /* REALIZAR EL CRUCE PARA Pm*Tp */
   for (i=0; i< size; i++)
     {
       if (!changed[i]) {
	  if (rate == 0) {
	     children[i] = mom[i];
	  }
	  else {
	     children[i] = dad[i];
	  }
	  continue;
       }

       m_domain->getValues(i, &min, &max);
       puntox=mom[i];  puntoy=dad[i];
       if (puntox > puntoy) {
	  swap(puntox,puntoy); 
       }

       I=(puntoy-puntox)*m_alpha;
       assert(I >= 0);

       A1 = puntox-I;

       if (A1 < min) 
	  A1 = min;

       B1 = puntoy+I;

	if (B1 > max) {
	    B1 = max;
	}

       children[i]=A1 + m_random->rand()*(B1-A1);
     }

}

 

CrossPBLX::CrossPBLX(double alpha) : m_alpha(alpha) {
}

CrossPBLX::~CrossPBLX(void) {
}

/**
 * Get the gene numbers to cross. 
 *
 * It grows lineally when current is increased, it starts from 0 when current == 0 to dim when current==maxtotal
 *
 * @param current current ratio of evaluations.
 * @param mintotal min ratio when it crosses all genes
 * @param dim gene dimension
 *
 * @return gene numbers to cross
 */
unsigned getNGen(double ratiocurrent, double mintotal, unsigned dim) {
   if (ratiocurrent >= mintotal) {
      return dim;
   }
   else {
      return (int)ceil(dim*(ratiocurrent/mintotal));
   }
}


void CrossPBLX::operator()(const tChromosomeReal &mom, tFitness fit_mom, const tChromosomeReal &dad, tFitness fit_dad, tChromosomeReal &children) {
    double puntox, puntoy, I, A1, B1;
    unsigned i, dim = dad.size();
    tGen min, max;

   assert(dad.size() == mom.size());

   // Obtain the gene numbers to cross 
   unsigned dim_ini = (unsigned) ceil(m_dim_ini*dim);
   unsigned dim_fin = (unsigned) ceil(m_dim_fin*dim);
   assert(dim_ini < dim_fin);

   // First copy all gen from mother
   copy(mom.begin(), mom.end(), children.begin());

   // Change the gene
   for (i=dim_ini; i< dim_fin; i++)
     {
       if (!m_domain->canBeChanged(i)) {
	    continue;
       }

       m_domain->getValues(i, &min, &max);
       puntox=mom[i];  puntoy=dad[i];
       I=fabs(puntoy-puntox)*m_alpha;

       A1 = mom[i]-I;

       if (A1 < min) 
	  A1 = min;

       B1 = mom[i]+I;

	if (B1 > max) {
	    B1 = max;
	}

       children[i]=A1 + m_random->rand()*(B1-A1);
     }

}
 
