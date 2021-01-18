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

#include "cmaeshan.h"
#include "origcmaes.h"
#include "cmaes_interface.h"
#include "distance.h"
#include "debug.h"
#include <algorithm>
#include <cassert>
#include <cstdio>

using namespace realea;
using namespace realea::internal;

typedef vector<tGen> DistVector;

class ProblemCMAESHansenEvalReal {
public:
    ProblemCMAESHansenEvalReal(IEval* eval) : m_eval(eval) {
    }

    tFitness eval(const double *sol, unsigned n) {
	 tChromosomeReal solv(n);

         for (unsigned i = 0; i < n; i++) {
	     solv[i] = sol[i];
         }
	 return m_eval->eval(solv);
   }

/*
    double eval(Real *sol, unsigned n) {
	 tChromosomeReal solv(n);

         for (unsigned i = 0; i < n; i++) {
	     solv[i] = sol[i];
         }
	 return m_eval->eval(solv);
   }
*/

private:
    IEval *m_eval;
};




class CMAESHansenParams : public ILSParameters {
public:
    cmaes_t evo;
    double *xinit, *stddev;
    double *fitvals;   /* objective function values of sampled population */
    bool init;

    CMAESHansenParams(unsigned dim) {
	xinit = new double[dim];
	stddev = new double[dim];
	fitvals = NULL;
	init = false;
    }

    ~CMAESHansenParams(void) {
	delete[] xinit;
	delete[] stddev;
	cmaes_exit(&evo);
    }
};

CMAESHansen::CMAESHansen(string fconfig) : m_debug(false) {
	m_rfactor = 0;
	m_nfactor = 0;
	m_fconfig = fconfig;
	m_bound_strategy = "always";
        m_setmyrandom = false;
        m_lambda = 0;
        m_mu = 0;
}

void CMAESHansen::setMyRandom(void) {
    m_setmyrandom = true;
}

/**
 * Obtain the range (max[0]-min[0], max[1]-min[1], ..., max[N-1]-min[N-1])
 *
 * @param domain
 *
 * @param range output
 */
static void getRange(DomainRealPtr domain, vector<tGen> *range) {
	tReal min, max;
	unsigned ndim = domain->getDimension();

	for (unsigned i = 0; i < ndim; i++) {
		domain->getValues(0, &min, &max);
		(*range)[i] = max-min;
	}
}


void CMAESHansen::searchRange(double factor) {
	assert(factor > 0 && factor <= 1);
	m_rfactor = factor;
}

void CMAESHansen::searchNeighborhood(double factor) {
	assert(factor > 0 && factor <= 1);
	m_nfactor = factor;
}

ILSParameters *CMAESHansen::getInitOptions(tChromosomeReal &sol) {
    int i, dim = sol.size();
    CMAESHansenParams *param;
    DistVector dist(dim);

   param = new CMAESHansenParams(dim);

   for (i = 0; i < dim; i++) {
      param->xinit[i] = sol[i];
   }

    if (m_nfactor) {
	if (m_pop == NULL) {
	   delete param;
	   throw ConfigException("CMAESHansen::Population");
	}

	min_vector_distance(sol, m_pop, dist);

	for (i = 0; i < dim; i++) {
	    param->stddev[i] = dist[i]*m_nfactor+0.001;
	}

   }
   else if (m_rfactor) {
	DomainRealPtr domain = m_problem->getDomain();
	vector<tGen> range(dim);
	getRange(domain, &range);

	for (i = 0; i < dim; i++) {
	    param->stddev[i] = range[i]*m_rfactor;
	}
   }

   return param;
}
    


void CMAESHansen::enableBoundChecking(void) {
   m_bound_strategy = "always";
}

void CMAESHansen::disableBoundChecking(void) {
   m_bound_strategy = "never";
}

void CMAESHansen::setPopsize(int lambda) {
    assert(lambda >= 0);
    m_lambda = lambda;
    
    if (m_lambda > 0 && m_lambda < m_mu) {
        m_mu = m_lambda;
        if (m_debug) {
           print_error("Warning: CMAES::mu is limited by CMAES::lambda to %d\n", m_mu);
        }
    }
}

void CMAESHansen::setParentsSize(int mu) {
    assert(mu >= 0);
    m_mu = mu;
    if (m_lambda > 0 && m_lambda < m_mu) {
        m_mu = m_lambda;
        if (m_debug) {
           print_error("Warning: CMAES::mu is limited by CMAES::lambda to %d\n", m_mu);
        }
    }
}

int CMAESHansen::getPopsize() {

    return m_lambda;
}

int CMAESHansen::getParentsSize() {

    return m_mu;
}

unsigned CMAESHansen::apply(ILSParameters *opt, tChromosomeReal &sol, tFitness &fitness, unsigned itera) {
   CMAESHansenParams *params = (CMAESHansenParams *) opt;
  double *const*pop; /* sampled population */
   unsigned dim;
  tFitness bestfit;
   int lambda = 0, counteval = 0;
  tFitness fbestever=0;
  double *xbestever=NULL; /* store best solution */
  double fmean;
  unsigned long seed;
  const double *xmean; 
   // User defined input parameters (need to be edited)
   char const * stop=NULL; /* stop message */
   int maxeval = itera-lambda;
   DomainRealPtr domain = m_problem->getDomain();

   if (m_setmyrandom) {
       seed = m_random->getSeed();
   }
   else {
       seed = 0;
   }

   //---------- BEGIN CMAESC ------------------------------------------
   ProblemCMAESHansenEvalReal *eval = new ProblemCMAESHansenEvalReal(m_eval);
   dim = sol.size();

   if (!params->init) {
      params->fitvals = cmaes_init(&params->evo, dim, params->xinit, params->stddev, seed, m_lambda, m_fconfig.c_str()); /* allocs fitvals */

      //if (m_lambda == 0) {
         m_lambda = cmaes_Get(&params->evo, "lambda");   /* needed for the restart */
      //}

      // cmaes_init sets params->evo.sp.mu to -1, so if we need to set this parameter we do it here, afterwards
      // if the given mu is within the bounds we use it, otherwise we use the autogenerated one
      if (m_mu > 0 && m_mu < m_lambda) {
         params->evo.sp.mu = m_mu;
      } else {
         m_mu = params->evo.sp.mu;
      }

      params->init = true;
   }

//   print_debug("cmaeshan.cc, params->evo.sp.mu: %d\n", params->evo.sp.mu);
//   print_debug("cmaeshan.cc, params->evo.sp.lambda: %d\n", params->evo.sp.lambda);

   lambda = m_lambda;
   maxeval = itera-lambda;

   bestfit = 0;
   params->evo.countevals = 0; /* a hack, effects the output and termination */

   while(!(stop=cmaes_TestForTermination(&params->evo)) && counteval < maxeval && !m_running->isFinish()) {
    	  /* Generate population of new candidate solutions */
	  pop = cmaes_SamplePopulation(&params->evo); /* do not change content of pop */

	  /* Here optionally handle constraints etc. on pop. You may
	   * call cmaes_ReSampleSingle(&params->evo, i) to resample the i-th
	   * vector pop[i], see below.  Do not change pop in any other
	   * way. You may also copy and modify (repair) pop[i] only
	   * for the evaluation of the fitness function and consider
	   * adding a penalty depending on the size of the
	   * modification.
	   */
	  /* Compute fitness value for each candidate solution */
	  for (int i = 0; i < cmaes_Get(&params->evo, "popsize") && !m_running->isFinish(); ++i) {
	    /* You may resample the solution i until it lies within the
	       feasible domain here, e.g. until it satisfies given  
               box constraints (variable boundaries). The function 
               is_feasible() needs to be user-defined.  
	       Assumptions: the feasible domain is convex, the optimum
	       is not on (or very close to) the domain boundary,
	       initialX is feasible and initialStandardDeviations are
	       sufficiently small to prevent quasi-infinite looping.
	    */
	    /* while (!is_feasible(pop[i])) 
	         cmaes_ReSampleSingle(&params->evo, i); 
	    */
	    // TODO: ¿poner?
	    //　cmaes_ReSampleSingle(&params->evo, i); 

	    if (m_bound_strategy == "always") {
		domain->clip(pop[i]);
	    }

	    params->fitvals[i] = eval->eval(pop[i], dim); 
	    counteval++;
            cmaes_print_info(&params->evo);

	    if (params->fitvals[i] < bestfit || counteval == 1) {
		bestfit = params->fitvals[i];
            }

	  }

	  /* update search distribution */
	  cmaes_UpdateDistribution(&params->evo, params->fitvals); 
	  
	  /* read control signals for output and termination */
//	  cmaes_ReadSignals(&params->evo, "signals.par"); /* from file signals.par */
	  
	  params->evo.countevals = counteval;     /* ditto */
	} /* while !cmaes_TestForTermination(&params->evo) */

	fbestever = cmaes_Get(&params->evo, "fbestever"); 
        xbestever = cmaes_GetInto(&params->evo, "xbestever", xbestever); /* alloc mem if needed */
      
	/* best estimator for the optimum is xmean, therefore check */
	xmean = cmaes_GetPtr(&params->evo, "xmean");

	if (counteval == 0) {
	    if (stop != NULL) {
		print_info("%s", stop);
	    }

	    if (fbestever) {
	    	copy(xbestever, xbestever+dim, sol.begin());
	    	fitness = fbestever;
	    }
            if (xbestever) {
                free(xbestever);
            }
	    delete eval;
	    return counteval;
	}


      if (!m_running->isFinish() && (fmean = eval->eval(xmean, dim)) < fbestever) {
	fbestever = fmean;
	xbestever = cmaes_GetInto(&params->evo, "xmean", xbestever);
      }

      if (fbestever < fitness) {
	copy(xbestever, xbestever+dim, sol.begin());
	fitness = fbestever;
      }
      if (xbestever) {
          free(xbestever);
      }
      delete eval;
      return counteval;
}

