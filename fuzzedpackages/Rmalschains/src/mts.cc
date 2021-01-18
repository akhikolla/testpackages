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
#include "mts.h"
#include "distance.h"
#include "random.h"
#include <cassert>
#include <cmath>

using namespace realea;

class MTSLSParams : public ILSParameters {
   public:
       MTSLSParams(unsigned gen_ini, double _delta) {
	  delta = _delta;
	  initialdelta = delta;
	  gen = gen_ini;
	  has_improved = true;
       }

       unsigned gen;
       bool has_improved;
       double delta;
       double initialdelta;

       ~MTSLSParams(void) {}

       virtual void store(tGen **aparams, unsigned *psize) {
	    unsigned size = 4;
	    tGen *params = new tGen[size];
	    params[0] = delta;
	    params[1] = initialdelta;
	    params[2] = gen;
	    params[3] = (has_improved ? 1.0 : 0.0);
	    *aparams = params;
	    *psize = size;
       }

       virtual void recover(tGen *params, unsigned size) {
	    assert(size > 3);
	    delta = params[0];
	    initialdelta = params[1];
	    gen = unsigned(params[2]);
	    
	    if (params[3] == 1) 
	       has_improved = true;
	    else
	       has_improved = false;
       }

};


MTSLS::MTSLS(double maxdelta, double mindelta, double ratio) {
    m_maxdelta = maxdelta;
    m_mindelta = mindelta;
    assert(ratio > 0 && ratio <= 1.0);
    m_ratio = ratio;
}

ILSParameters *MTSLS::recoverOptions(tGen *params, unsigned size) {
    MTSLSParams *option;
    DomainRealPtr domain = m_problem->getDomain();
    unsigned ndim = domain->getDimension();
    unsigned gen_ini;
 
    for (gen_ini=0; gen_ini < ndim && !domain->canBeChanged(gen_ini); ++gen_ini) 
	;

    assert(gen_ini < ndim);

    option = new MTSLSParams(gen_ini, 0);
    option->recover(params, size);
    return ((ILSParameters *) option);
}

void MTSLS::storeOptions(ILSParameters *params, tGen **paparams, unsigned *psize) {
     unsigned size = 1;

     if (params != NULL) {
	MTSLSParams *p = (MTSLSParams *) params;
	p->store(paparams, psize);
	assert(size == *psize);
     }
     else {
	*paparams = NULL;
     }

     *psize = size;
}


ILSParameters *MTSLS::getInitOptions(tChromosomeReal &sol) {
    DomainRealPtr domain = m_problem->getDomain();
    unsigned gen_ini, ndim;

    ndim = sol.size();
    MTSLSParams *option;
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

    for (gen_ini=0; gen_ini < ndim && !domain->canBeChanged(gen_ini); ++gen_ini) 
	;

    assert(gen_ini < ndim);

    option = new MTSLSParams(gen_ini, step);
    return ((ILSParameters *) option);
}

unsigned MTSLS::apply(ILSParameters *params, tChromosomeReal &sol, tFitness &sol_perf, unsigned maxeval) {
   MTSLSParams *p = (MTSLSParams *) params;
   unsigned ndim = sol.size();
   DomainRealPtr domain = m_problem->getDomain();
   tFitness newsol_perf, sign;
   tGen oldgen;
   unsigned numEval = 0;

   for (numEval = 0; numEval < maxeval && !m_running->isFinish(); ) {

      if (p->gen == 0 && !p->has_improved) {
	 p->delta /= 2;

	 if (p->delta < m_mindelta) {
	    p->delta = p->initialdelta;
	 }
      }

      while(p->gen < ndim && numEval < maxeval && !m_running->isFinish() ) {
	  sign = (m_random->rand() > 0.5 ? 1 : -1);

	  if (m_random->rand() > m_ratio) {
	     continue;
	  }

	  oldgen = sol[p->gen];
	  sol[p->gen] = domain->clip(p->gen, oldgen+sign*p->delta);
	  numEval++;
	  newsol_perf = m_eval->eval(sol);

	// Si lo mejoro lo copio
	if (m_problem->isBetter(newsol_perf, sol_perf) ) {
	    sol_perf = newsol_perf;
	    p->has_improved = true;
	}
	else { 
	    sol[p->gen] = oldgen;

	   if (numEval < maxeval && !m_problem->isBetter(newsol_perf, sol_perf) && !m_running->isFinish()) {
	    sol[p->gen] = domain->clip(p->gen, oldgen-sign*0.5*p->delta);
	    newsol_perf = m_eval->eval(sol);
	    numEval++;
	    
	    if (m_problem->isBetter(newsol_perf, sol_perf)) {
		sol_perf = newsol_perf;
		p->has_improved = true;
	    }
	    else  
		sol[p->gen] = oldgen;
	    } 

	}

	do {
	    p->gen = (p->gen + 1) % ndim;
	
	    if (p->gen == 0) {
		p->has_improved = false;
	    }

	} while (!domain->canBeChanged(p->gen));


      } // Iterate all dimension


   } // De comprobar las iteraciones

   return numEval;
}
