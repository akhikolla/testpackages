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
#include "mts2.h"
#include "distance.h"
#include "random.h"
#include <cassert>
#include <cmath>

using namespace realea;

class MTSLS2Params : public ILSParameters {
   public:
       MTSLS2Params(double _delta) {
	  delta = _delta;
	  initialdelta = delta;
	  has_improved = true;
       }

       bool has_improved;
       double delta;
       double initialdelta;

       ~MTSLS2Params(void) {}

       virtual void store(tGen **aparams, unsigned *psize) {
	    unsigned size = 3;
	    tGen *params = new tGen[size];
	    params[0] = delta;
	    params[1] = initialdelta;
	    params[2] = (has_improved ? 1.0 : 0.0);
	    *aparams = params;
	    *psize = size;
       }

       virtual void recover(tGen *params, unsigned size) {
	    assert(size > 3);
	    delta = params[0];
	    initialdelta = params[1];
	    
	    if (params[2] == 1) 
	       has_improved = true;
	    else
	       has_improved = false;
       }

};


MTSLS2::MTSLS2(double maxdelta, double mindelta) {
    m_maxdelta = maxdelta;
    m_mindelta = mindelta;
}

ILSParameters *MTSLS2::recoverOptions(tGen *params, unsigned size) {
    MTSLS2Params *option;

    option = new MTSLS2Params(0);
    option->recover(params, size);
    return ((ILSParameters *) option);
}

void MTSLS2::storeOptions(ILSParameters *params, tGen **paparams, unsigned *psize) {
     unsigned size = 1;

     if (params != NULL) {
	MTSLS2Params *p = (MTSLS2Params *) params;
	p->store(paparams, psize);
	assert(size == *psize);
     }
     else {
	*paparams = NULL;
     }

     *psize = size;
}


ILSParameters *MTSLS2::getInitOptions(tChromosomeReal &sol) {
    MTSLS2Params *option;
    unsigned nearest;
    /**
     * Calculo sigma como la mitad de la distancia al mas cercano
     */
    double dist;

    if (m_pop == NULL) 
       dist = 0.2;
    else 
       dist = distanceMin(sol, m_pop, &nearest);

    double step = dist/2.0;

    if (step > m_maxdelta) {
	step = m_maxdelta;
    }

    option = new MTSLS2Params(step);
    return ((ILSParameters *) option);
}

unsigned MTSLS2::apply(ILSParameters *params, tChromosomeReal &sol, tFitness &sol_perf, unsigned maxeval) {
   MTSLS2Params *p = (MTSLS2Params *) params;
   unsigned ndim = sol.size();
   tChromosomeReal oldsol(ndim);
   unsigned *D = new unsigned[ndim];
   unsigned *r = new unsigned[ndim];
   DomainRealPtr domain = m_problem->getDomain();
   tFitness newsol_perf;
   unsigned numEval = 0;
   unsigned d;

   if (!p->has_improved) {
      p->delta /= 2;

      if (p->delta < m_mindelta) {
	 p->delta = p->initialdelta;
      }
   }

   p->has_improved = false;

   while(numEval < maxeval && !m_running->isFinish() ) {
      oldsol = sol;

      for (d = 0; d < ndim; d++) {
	 D[d] = -1 + 2*m_random->randint(0, 1);
	 r[d] = m_random->randint(0, 3);
      }

      for (d = 0; d < ndim; d++) {
	 if (r[d] == 0) {
	    sol[d] = domain->clip(d, sol[d]-D[d]*p->delta);
	 }
      }

      numEval++;
      newsol_perf = m_eval->eval(sol);

      // Si lo mejoro lo copio
      if (m_problem->isBetter(newsol_perf, sol_perf) ) {
	 sol_perf = newsol_perf;
	 p->has_improved = true;
      }
      else { 
	 sol = oldsol;

	 if (numEval < maxeval && !m_running->isFinish()) {
	    for (d = 0; d < ndim; d++) {
	       if (r[d] == 0) {
		  sol[d] = domain->clip(d, sol[d]+0.5*D[d]*p->delta);
	       }
	    }

	    newsol_perf = m_eval->eval(sol);
	    numEval++;

	    if (m_problem->isBetter(newsol_perf, sol_perf)) {
	       sol_perf = newsol_perf;
	       p->has_improved = true;
	    }
	    else 
	      sol = oldsol;
	 } 

      }

   } // Iterate all dimension

   delete[] D;
   delete[] r;
   return numEval;
}
