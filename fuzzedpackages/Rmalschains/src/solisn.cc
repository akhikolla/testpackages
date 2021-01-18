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

#include "solisn.h"
#include "distance.h"
#include "random.h"
#include <functional>
#include <algorithm>
#include <cassert>
#include <cmath>

using namespace realea;

class SWDimParams : public ILSParameters {
   public:
       SWDimParams(unsigned dim) : bias(dim) {
       }
       unsigned dim(void) {
	  return bias.size();
       }

       double delta;
       vector<double> bias;
       unsigned numFailed;
       unsigned numSuccess;
       virtual ~SWDimParams(void) {}

       virtual void store(tGen **aparams, unsigned *psize) {
	    unsigned size = 3+bias.size();
	    tGen *params = new tGen[size];
	    params[0] = delta;
	    params[1] = numFailed;
	    params[2] = numSuccess;
	    copy(bias.begin(), bias.end(), params+3);
	    *aparams = params;
	    *psize = size;
       }

       virtual void recover(tGen *params, unsigned size) {
	    assert(size > 3);
	    delta = params[0];
	    numFailed = (unsigned) ceil(params[1]);
	    numSuccess = (unsigned) ceil(params[2]);
	    copy(params+3, params+size, bias.begin());
       }

};


SWNDim::SWNDim(void) {
    m_maxdelta = m_mindelta = -1;
}

void SWNDim::setDelta(double maxdelta) {
   assert(maxdelta > 0); 
   m_maxdelta = maxdelta;
}

void SWNDim::setDelta(double mindelta, double maxdelta) {
    assert(mindelta > 0 && mindelta <= maxdelta);
    m_mindelta = mindelta;
    m_maxdelta = maxdelta;
}

ILSParameters *SWNDim::getInitOptions(tChromosomeReal &sol) {
    SWDimParams *option;
    unsigned dim = m_problem->getDimension();
    option = new SWDimParams(dim);
    option->numFailed = option->numSuccess = 0;
    unsigned nearest;
    if (m_pop == NULL) {
       assert(m_maxdelta > 0);
       option->delta = m_maxdelta;
    }
    else {
       /**
	* Calculo sigma como la mitad de la distancia al mas cercano
	*/
       double dist = distanceMin(sol, m_pop, &nearest);
       option->delta = dist/2.0;

       if (m_maxdelta > 0 && option->delta > m_maxdelta) {
	  option->delta = m_maxdelta;
       }
       else if (m_mindelta > 0 && option->delta < m_mindelta) {
	  option->delta = m_mindelta;
       }

    }

    fill(option->bias.begin(), option->bias.end(), 0.0);
    
    return ((ILSParameters *) option);
}

tFitness SWNDim::getNeighbour(ILSParameters *param, tChromosomeReal &actual, tChromosomeReal &dif, tChromosomeReal &newsol, bool *change) {
    unsigned i;
    SWDimParams *p = (SWDimParams *) param;
    unsigned ndim = actual.size();
    DomainRealPtr domain = m_problem->getDomain();

    for (i = 0; i < ndim; i++) {
        if (change[i]) {
	    dif[i] = m_random->normal(p->delta);
	}
	else {
	    dif[i] = 0;
	}

	newsol[i] = actual[i] + p->bias[i] +dif[i];
    }

    domain->clip(newsol);
    return m_eval->eval(newsol);
}



static double increm_bias(const double &bias, const double &dif) {
       return 0.2*bias+0.4*(dif+bias);
}

static double dec_bias(const double &bias, const double &dif) {
    	return bias-0.4*(dif+bias);
}

static double divide_bias(const double &bias, const double &dif) {
    	return bias*0.5;
}

static unsigned maxChanged(unsigned m_strategy, unsigned ndim) {
    if (m_strategy == 0)
      return 0;
    else if (m_strategy == 1) {
      return unsigned(0.2*ndim <= 10 ? 0.2*ndim : 10);
    }
    else if (m_strategy == 5) {
      return unsigned(0.2*ndim <= 05 ? 0.2*ndim : 50);
    }

    return 0;
}


static void setChanged(unsigned m_strategy, Random *m_random, bool *changed, unsigned ndim) {
   int n = m_random->randint(0, ndim-1);
   unsigned int i;
/*
   for (i = 0; i < ndim && changed[n]; i++) {
	n = (n+1)%ndim;
   }
*/

   fill(changed, changed+ndim, false);

   if (m_strategy == 0) {
   }
   else if (m_strategy == 1 || m_strategy == 5) {
     unsigned min = maxChanged(m_strategy, ndim);

     for (i = 0; i < min; i++) {
	changed[(n+i) %ndim] = true;
     }

   }
   else if (m_strategy == 2) {
   	for (i = 0; i < ndim && m_random->rand() < 0.9; i++) {
		changed[(n+i) %ndim] = true;
	}

   }
   else if (m_strategy == 3) {
   	for (i = 0; i < ndim; i++) {
	    if (m_random->rand() < 0.1) {
		changed[(n+i) %ndim] = true;
	    }
	}

   }

}



void SWNDim::setStrategy(unsigned st) {
	m_strategy = st;
}

unsigned SWNDim::apply(ILSParameters *params, tChromosomeReal &sol, tFitness &sol_perf, unsigned maxeval) {
   SWDimParams *p = (SWDimParams *) params;
   unsigned ndim = sol.size();
   tChromosomeReal dif(ndim), newsol(ndim);
   DomainRealPtr domain = m_problem->getDomain();
   bool *changed = new bool[ndim];
   tFitness newsol_perf;
   unsigned gen, itera;
   unsigned numEval = 0;

   if (m_strategy <= 3) {
   	setChanged(m_strategy, m_random,changed, ndim); 
   }
   if (m_strategy == 6) {
   	setChanged(1, m_random,changed, ndim); 
   }
   
   itera = 1;
   unsigned iteradim = 100;
   unsigned maxitera7 = maxeval/10;
   unsigned maxitera8 = maxeval/5;

   for (numEval = 0; numEval < maxeval && !m_running->isFinish(); itera++) {
      if (m_strategy == 5 && ( (itera % iteradim)==0)) {
	 setChanged(2, m_random, changed, ndim);
      }
      else if (m_strategy == 6 && itera > 0 && ( (itera % iteradim/2)==0)) {
	 setChanged(1, m_random, changed, ndim);
      }
      else if (m_strategy == 7 && itera > 0 && ( (itera % maxitera7/2)==0)) {
	 setChanged(5, m_random, changed, ndim);
      }
      else if (m_strategy == 8 && itera > 0 && ( (itera % maxitera8/2)==0)) {
	 setChanged(5, m_random, changed, ndim);
      }



      newsol_perf = getNeighbour(p, sol, dif, newsol, changed);
      numEval++;

      // Si lo mejoro lo copio
      if (m_problem->isBetter(newsol_perf, sol_perf) ) {
	 copy(newsol.begin(), newsol.end(), sol.begin());
	 sol_perf = newsol_perf;

	 // Adapto bias 
	 transform(p->bias.begin(), p->bias.end(), dif.begin(), p->bias.begin(), 
		     increm_bias);
	 p->numSuccess++;
	 p->numFailed = 0;
      }
      else if (numEval < maxeval && !m_problem->isBetter(newsol_perf, sol_perf) && !m_running->isFinish()) {
	 // Invierto la dirección
	 for (gen = 0; gen < ndim; gen++) {
	     newsol[gen] = sol[gen] - p->bias[gen] - dif[gen];
	 }

	 domain->clip(newsol);
	 newsol_perf = m_eval->eval(newsol);
	 numEval++;

	 if (m_problem->isBetter(newsol_perf, sol_perf)) {
	    copy(newsol.begin(), newsol.end(), sol.begin());
	    sol_perf = newsol_perf;

	    // Disminuyo bias
	    transform(p->bias.begin(), p->bias.end(), dif.begin(), p->bias.begin(), 
			dec_bias);

	    p->numSuccess++;
	    p->numFailed = 0;
	 }
	 else {
	    // Se reduce bias a la mitad
	    transform(p->bias.begin(), p->bias.end(), dif.begin(), p->bias.begin(), 
			divide_bias);
	    p->numFailed++;
	    p->numSuccess = 0;
	 }

      }

      if (p->numSuccess >= maxSuccess) {
	 p->numSuccess = 0;
	 p->delta *= 2;
      }
      else if (p->numFailed >= maxFailed) {
	 p->numFailed = 0;
	 p->delta *= 0.5;
      }

   } // De comprobar las iteraciones

   delete[] changed;
//   printf("p->delta: %e -> %e\n", oldelta, p->delta);
   return numEval;
}
