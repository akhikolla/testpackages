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

#include "solisn2.h"
#include "distance.h"
#include "random.h"
#include <functional>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>

using namespace realea;

class SW2DimParams : public ILSParameters {
   public:
       SW2DimParams(unsigned dim) : delta(dim), delta_init(dim), bias(dim) {
       }
       unsigned dim(void) {
	  return bias.size();
       }

       vector<double> delta;
       vector<double> delta_init;
       vector<double> bias;
       unsigned numFailed;
       unsigned numSuccess;
       virtual ~SW2DimParams(void) {}
};


SWN2Dim::SWN2Dim(void) {
    m_maxdelta = m_mindelta =  m_initialDelta = -1;
}

void SWN2Dim::setDelta(double maxdelta) {
   assert(maxdelta > 0); 
   m_maxdelta = maxdelta;
}

void SWN2Dim::seInitialtDelta(double init_delta){
    m_initialDelta = init_delta;
}

void SWN2Dim::setDelta(double mindelta, double maxdelta) {
    assert(mindelta > 0 && mindelta <= maxdelta);
    m_mindelta = mindelta;
    m_maxdelta = maxdelta;
}

ILSParameters *SWN2Dim::getInitOptions(tChromosomeReal &sol) {
    SW2DimParams *option;
    unsigned dim = m_problem->getDimension();
    option = new SW2DimParams(dim);
    option->numFailed = option->numSuccess = 0;
    unsigned nearest;

    if(m_initialDelta > 0){
        //printf("delta according to distance to region\n");
        fill(option->delta.begin(), option->delta.end(), m_initialDelta);
        fill(option->delta_init.begin(), option->delta_init.end(), m_initialDelta);
    }
    else if (m_pop == NULL) {
       assert(m_maxdelta > 0);

       //printf("delta with max delta \n");
       fill(option->delta.begin(), option->delta.end(), m_maxdelta);
       fill(option->delta_init.begin(), option->delta_init.end(), m_maxdelta);
    }
    else {
       /**
	* Calculo sigma como la mitad de la distancia al mas cercano
	*/
       //printf("delta according to distance to neighbour\n");
       for (unsigned i = 0; i < dim; i++) {
	    distanceMin(sol, m_pop, &nearest);
	    double dist = fabs(sol[i] - m_pop->getInd(nearest)->gen(i));
	    option->delta_init[i] = dist/2.0;

	    if (m_maxdelta > 0 && option->delta[i] > m_maxdelta) {
		option->delta_init[i] = m_maxdelta;
	    }
	    else if (m_mindelta > 0 && option->delta[i] < m_mindelta) {
		option->delta_init[i] = m_mindelta;
	    }

       }

       copy(option->delta_init.begin(), option->delta_init.end(), option->delta.begin());

    }

    fill(option->bias.begin(), option->bias.end(), 0.0);

    return ((ILSParameters *) option);
}

tFitness SWN2Dim::getNeighbour(ILSParameters *param, tChromosomeReal &actual, tChromosomeReal &dif, tChromosomeReal &newsol, bool *change) {
   unsigned i;
   SW2DimParams *p = (SW2DimParams *) param;
   unsigned ndim = actual.size();
   DomainRealPtr domain = m_problem->getDomain();

   for (i = 0; i < ndim; i++) {
      if (change[i]) {
	 dif[i] = m_random->normal(p->delta[i]);
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
   if (dif == 0)
     return bias;
   else
   return bias-0.4*(dif+bias);
}

static double divide_bias(const double &bias, const double &dif) {
   if (dif == 0)
     return bias;
   else
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
      unsigned int min = maxChanged(m_strategy, ndim);

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



void SWN2Dim::setStrategy(unsigned st) {
   m_strategy = st;
}

unsigned SWN2Dim::apply(ILSParameters *params, tChromosomeReal &sol, tFitness &sol_perf, unsigned maxeval) {
   SW2DimParams *p = (SW2DimParams *) params;
   unsigned ndim = sol.size();
   tChromosomeReal dif(ndim), newsol(ndim);
   DomainRealPtr domain = m_problem->getDomain();
   bool *changed = new bool[ndim];
   tFitness newsol_perf;
   unsigned gen, itera, i;
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

	 for (i = 0; i < ndim; ++i) {
	    if (changed[i])	
		p->delta[i] *= 2;

            if (p->delta[i] > m_maxdelta){
                //printf("max delta reached:\n");
		p->delta[i] = m_maxdelta;
            }
	 }
      }
      else if (p->numFailed >= maxFailed) {
	 p->numFailed = 0;

	 for (i = 0; i < ndim; ++i) {
	    if (changed[i]) {
		p->delta[i] *= 0.5;

                if (p->delta[i] < m_mindelta){
                    //printf("min delta reached:\n");
		  p->delta[i] = p->delta_init[i];
                }
	    }
	 }
      }

   } // De comprobar las iteraciones
   delete[] changed;
   return numEval;
}
