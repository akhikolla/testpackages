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

#include "cmaes.h"
#include "cmaesb.h"
#include "distance.h"
#include <algorithm>
#include <cassert>

#define _PRINTDEBUG 0

using namespace realea;
using namespace realea::internal;

class ProblemEvalReal : public IEvalReal {
public:
    ProblemEvalReal(IEval* eval) : m_eval(eval) {
    }

    double eval(Real *sol, unsigned n) {
	 tChromosomeReal solv(n);

         for (unsigned i = 0; i < n; i++) {
	     solv[i] = sol[i];
         }
	 return m_eval->eval(solv);
   }

private:
    IEval *m_eval;
};

void copySol(const Real *ini, const Real *fin, tChromosomeReal &sol) {
   const Real *p;
   unsigned i;

   for (i = 0, p = ini; p != fin; ++p, ++i) {
	sol[i] = *p;
   }

}

typedef vector<tGen> DistVector;

class CMAESParams : public ILSParameters {
public:
    MyMatrix B;
    MyMatrix BD;
    // covarianza matrix
    MyMatrix C;
    DiagonalMatrix D;
    int m_N;

    // Initialize dynamic (internal) strategy parameters and constants
    ColumnVector pc;
    // evolution paths for C and sigma
    ColumnVector ps;
    // coordinate wise standard deviation (step size)
    ColumnVector sigma;

    CMAESParams(int N, ColumnVector &psigma) : B(N,N), BD(N,N), C(N,N), D(N), pc(N), ps(N), sigma(N) {
       m_N = N;
       reset(psigma);
    }

    void reset(ColumnVector &sigma) {
        B = eye(m_N);
	D << eye(m_N);
	BD = (B*D);
	C = (BD* BD.t()); // covariance matrix
	pc = 0;
	ps = 0;
	this->sigma = sigma;
    }

};

MyReturnMatrix randn(Random *m_random, int fil,int col) {
    int total = fil*col;
    double *data = new double[total];
    MyMatrix R(fil,col);
    int i;

    for (i = 0; i < total; i++) {
	data[i] = m_random->normal(1);
    }

    R <<data;
    R.Release();
    delete[] data;
    return R;
}



CMAES::CMAES(void) : m_debug(false) {
	m_pop = NULL;
	m_rfactor = 0;
	m_nfactor = 0;
}

/**
 * Obtain the range (max[0]-min[0], max[1]-min[1], ..., max[N-1]-min[N-1])
 *
 * @param domain
 *
 * @param range output
 */
void getRange(DomainRealPtr domain, vector<tGen> *range) {
	tReal min, max;
	unsigned ndim = domain->getDimension();

	for (unsigned i = 0; i < ndim; i++) {
		domain->getValues(0, &min, &max);
		(*range)[i] = max-min;
	}
}


double getMean(DomainRealPtr domain) {
	double sum=0;
	tReal min, max;
	unsigned ndim = domain->getDimension();

	for (unsigned i = 0; i < ndim; i++) {
		domain->getValues(0, &min, &max);
		sum += (max-min);
	}

	return sum/ndim;
}

void CMAES::searchRange(double factor) {
	assert(factor > 0 && factor <= 1);
	m_rfactor = factor;
}

void CMAES::searchNeighborhood(double factor) {
	assert(factor > 0 && factor <= 1);
	m_rfactor = factor;
}


ILSParameters *CMAES::getInitOptions(tChromosomeReal &sol) {
   unsigned N = sol.size();
   DistVector dist(N);
   ColumnVector sigma(N);

   if (m_nfactor) {
	if (m_pop == NULL) {
	   throw ConfigException("CMAES::Population");
	}

	min_vector_distance(sol, m_pop, dist);
	copyToColumn(dist, &sigma);
	sigma *= m_nfactor;
   }
   else if (m_rfactor) {
	DomainRealPtr domain = m_problem->getDomain();
	vector<tGen> range(N);
	getRange(domain, &range);
	copyToColumn(range, &sigma);
	sigma *= m_rfactor;
   }
   return new CMAESParams(m_problem->getDimension(), sigma);
}

unsigned CMAES::apply(ILSParameters *opt, tChromosomeReal &sol, double &fitness, unsigned itera) {
   ColumnVector better;
   CMAESParams *params = (CMAESParams *) opt;
   double fit_better, fit_actual;
   // User defined input parameters (need to be edited)
   int N = m_problem->getDimension();
   unsigned maxeval = itera;
   DomainRealPtr domain = m_problem->getDomain();

   //---------- BEGIN CMAESC ------------------------------------------
   IEvalReal *eval = new ProblemEvalReal(m_eval);

   ////  Init in avanced mode 
   CMAESBound bound(eval, domain);

   // Active the bound checks if the bounds are valids
   ColumnVector xmean(N);

   copyToColumn(sol, &xmean);

   better = xmean;
   fit_better = fitness;

   int stopeval = maxeval;         // stop after stopeval number of function evaluations

   // Strategy parameter setting: Selection  
   int lambda=4+(int)floor(3*log(N));  // population size, offspring number
   int mu = (int) floor(lambda/2.0);      // number of parents/points for recombination

   // Indices para escoger los n mejores
   int *arindex = new int[lambda];
   int *arindex_raw = new int[lambda];

   // muXone array for weighted recombination
   ColumnVector weights(mu);
   int i;

   for (i = 1; i <= mu; i++) {
      weights[i-1] = log(mu+1) - log(i);
   }

   double mueff=pow2(weights.Sum())/pow2(weights).Sum(); //variance-effective size of mu

#if _PRINTDEBUG
   if (m_debug) {
      cout <<"xmean: \n" <<xmean <<endl;
      cout <<"C: \n" <<params->C <<endl;
      cout <<"B: \n" <<params->B <<endl;
      cout <<"sigma: " <<params->sigma <<endl;
      cout <<"Lambda: " <<lambda <<endl;
      cout <<"Mu: " <<mu <<endl;
      cout <<"Weights: \n" <<weights <<endl;
      cout <<"mueff " <<mueff <<endl;
   }
#endif

   // Strategy parameter setting: Adaptation
   double cc = 4.0/(N+4.0);    // time constant for cumulation for covariance matrix
   double cs = (mueff+2)/(N+mueff+3);   // t-const for cumulation for sigma control
   double mucov = mueff;   // size of mu used for calculating learning rate ccov
   double ccov = (1/mucov) * 2/pow2(N+1.4) + (1-1/mucov)      // learning rate for
     *((2*mueff-1)/(pow2(N+2)+2*mueff));            // covariance matrix

   double damps = 1 + 2*max(0.0, sqrt((mueff-1)/(N+1))-1) + cs; // damping for sigma 
   // usually close to 1
#if _PRINTDEBUG
   if (m_debug) {    
      cout <<"cc: " <<cc <<endl;
      cout <<"cs: " <<cs <<endl;
      cout <<"ccov: " <<ccov <<endl;
      cout <<"damps: " <<damps <<endl;
   }
#endif 

   DiagonalMatrix DiagZeros(N);
   DiagZeros = 0;

   double chiN=sqrt(N)*(1-1/(4.0*N)+1/(21.0*pow2(N)));  // expectation of 
   //   ||N(0,I)|| == norm(randn(N,1))
   //

   weights = weights/weights.Sum();     // normalize recombination weights array

#if _PRINTDEBUG
   if (m_debug) {
      cout <<" C : \n" <<params->C <<endl;
      cout <<"chiN : " <<chiN <<endl;
      cout <<"weights : \n" <<weights <<endl;
   }
#endif

   // -------------------- Generation Loop --------------------------------

   int counteval = 0;  
   RowVector arfitness(lambda);
   RowVector arfitness_sel(lambda);
   RowVector arfitness_raw(lambda);
   MyMatrix arx(N,lambda);
   MyMatrix arxvalid(N,lambda);
   Matrix arx_bests(N,mu);
   Matrix arz_bests(N,mu);

   // Check if the current step is the last
   bool last = false;
   int countiter = 0;

   while (counteval+lambda < stopeval && !m_running->isFinish()  && !last 
	       // && !converged
	 ) {
      countiter += 1;
      // Generate and evaluate lambda offspring
      MyMatrix arz = randn(m_random,N,lambda);  // array of normally distributed mutation vectors

#if _PRINTDEBUG
      if (m_debug) {	
	 cout <<"randn:\n" <<arz <<endl;
	 cout <<"xmean:\n" <<xmean <<endl;
	 cout <<"sigma:\n" <<params->sigma <<endl;
	 cout <<"BD:\n" <<params->BD <<endl;
      }
#endif

      int k;
 
      // Set
      bound.setParam(lambda, mueff, params->sigma, params->C);

      for (k = 1; k <= lambda && counteval < stopeval; k++) {
	 ColumnVector arx_k = xmean+DotVectors(params->sigma, (params->BD*arz.Column(k))); // add mutation Eq. (1)
	 arx.Column(k) <<arx_k;
	 counteval += 1;
      }

      // Check k >= lambda 
      if (k < lambda) {
	 lambda = min(k,lambda);
	 last = true;
      }

      // Evaluo los nuevos puntos, obtengo dos vectores de fitness
      // el primero contiene el fitness real y el segundo el empleado para
      // la selección de la media
      bound.evalSols(xmean, arx, arxvalid, arfitness_raw, arfitness_sel);

      arfitness = arfitness_sel;

      if (counteval == stopeval) {
	 last = true;
      }

#if _PRINTDEBUG
      if (m_debug) {	
	 cout <<"arz :\n" <<arz <<endl;
	 cout <<"arx :\n" <<arx <<endl;
	 cout <<"arxvalid :\n" <<arxvalid <<endl;
	 cout <<" arfitness :\n" <<arfitness <<endl;
      }
#endif

      range(1, lambda, arindex);
      set_sort_matrix(&arfitness);
      partial_sort(arindex, arindex+mu, arindex+lambda, sort_index_matrix);

      range(1, lambda, arindex_raw);
      set_sort_matrix(&arfitness_raw);
      partial_sort(arindex_raw, arindex_raw+1, arindex_raw+lambda, sort_index_matrix);


      // ReActualizo el mejor obtenido hasta el momento
      int newbest = arindex_raw[0];

#if _PRINTDEBUG
      if (m_debug) {	
	 cout <<"arfitness[newbest]" <<arfitness_raw[newbest-1] <<endl;
      }
#endif

      fit_actual = arfitness_raw[newbest-1];

      if (m_running->isBetter(fit_actual, fit_better)) {
	 fit_better = fit_actual;
	 better = arxvalid.Column(newbest);
      }

      if (last) {
	 copySol(better.Store(), better.Store()+N, sol);
	 fitness = fit_better;
	 return counteval;
      }

      // Compruebo la convergencia
      int posconv = (int) ceil(1.0+lambda/4.0);

      if (arfitness[arindex[0]-1] == arfitness[arindex[posconv]-1]) {
	 //	   cout <<"CMAESC: Converged at " <<counteval <<endl;
	 //	   converged = true;
	 params->sigma *= exp(0.2+cs/damps);
      }

      // Update the local search parameters

      //Sort by fitness and compute weighted mean into xmean
      getColumns(arxvalid, arindex, mu, arx_bests);
      getColumns(arz, arindex, mu, arz_bests);

      xmean = arx_bests*weights;   // recombination, new mean value
      MyMatrix zmean = arz_bests*weights;   // == sigma^-1*D^-1*B'*(xmean-xold)
      //  Cumulation: Update evolution paths
      params->ps = (1-cs)*params->ps + sqrt(cs*(2-cs)*mueff) * (params->B * zmean);            // Eq. (4)

      double hsig = norm(params->ps)/pow(sqrt(1-(1-cs)), ((2.0*counteval)/lambda))/chiN < 1.5 + 1/(N+1);

      params->pc = (1-cc)*params->pc + hsig * sqrt(cc*(2-cc)*mueff) * (params->B * params->D * zmean); // Eq. (2)

#if _PRINTDEBUG
      if (m_debug) {	
	 cout <<"arx_bests\n" <<arx_bests <<endl;
	 cout <<"arz_bests\n" <<arz_bests <<endl;
	 cout <<"new xmean:\n" <<xmean <<endl;
	 cout <<"zmean:\n" <<zmean <<endl;
	 cout <<"ps:\n" <<params->ps <<endl;
	 cout <<"hsig:\n" <<hsig <<endl;
	 cout <<"pc:\n" <<params->pc <<endl; 
      }
#endif

      // Adapt covariance matrix C
      params->C = (1-ccov) * params->C                     // regard old matrix      // Eq. (3)
	+ ccov * (1/mucov) * (params->pc*params->pc.t()     // plus rank one update
		    + (1-hsig) * cc*(2-cc) * params->C) 
	+ ccov * (1-1/mucov) //...           // plus rank mu update 
	* (params->B*params->D*arz_bests) // ...
	*  weights.AsDiagonal() * (params->BD*arz_bests).t();               

#if _PRINTDEBUG
      if (m_debug) {
	 cout <<"C:\n" <<params->C <<endl;
      }
#endif

      // Adapt step size sigma
      params->sigma = params->sigma * exp((cs/damps)*(norm(params->ps)/chiN - 1));             // Eq. (5)


//       if (ls->checkLocal() && params->sigma > maxsigma) {
	 //	    cout <<"CMAESC Warning: Sigma value too high at " <<counteval <<endl;
//	 last = true;
//      }

      // Update B and D from C
      // This is O(N^3). When strategy internal CPU-time is critical, the
      // next three lines can be executed only every (alpha/ccov/N)-th
      // iteration step, where alpha is e.g. between 0.1 and 10 

      UpperTriangularMatrix triu,triu1;

      triu << params->C; triu1 <<params->C;
      triu1.Inject(DiagZeros);

#if _PRINTDEBUG
      if (m_debug) {	
	 cout <<"triu:\n" <<triu <<endl;
	 cout <<"triu1:\n" <<triu1 <<endl;
      }
#endif

      params->C = triu+triu1.t(); // Enforce symmetric

      SymmetricMatrix S(N);
      S <<params->C; 
#if _PRINTDEBUG
      if (m_debug) {	
	 cout <<"C:\n" <<params->C <<endl;
	 cout <<"S:\n" <<S <<endl;
      }
#endif

      EigenValues(S,params->D,params->B);  // eigen decomposition, B==normalized eigenvectors

#if _PRINTDEBUG
      if (m_debug) {	
	 cout <<"sigma:\n" <<params->sigma <<endl;
	 cout <<"C:\n" <<params->C <<endl;
	 cout <<"B: \n" <<params->B <<endl;
	 cout <<"D: \n" <<params->D <<endl;
      }
#endif

      // Compruebo que explore en toda dimensión
      checkDiag(params->C, params->D);

      params->D << sqrt(params->D);  // D contains standard deviations now

#if _PRINTDEBUG
      if (m_debug) {
	 cout <<"Sqrt D: \n" <<params->D <<endl;
	 cout <<counteval <<": " <<arfitness[newbest-1] <<endl;
      }
#endif
      // Its value is used again
      params->BD = (params->B*params->D);

      // Compruebo los ejes
      checkAxis(xmean, ccov, cs, damps, countiter, params->sigma, params->C, params->BD);
   }

#if _PRINTDEBUG
   if (m_debug) {
      cout <<"xmean: \n" <<xmean <<endl;
      cout <<"better: \n" <<better <<endl;
      cout <<"C: \n" <<params->C <<endl;
      cout <<"B: \n" <<params->B <<endl;
      cout <<"sigma: " <<params->sigma <<endl;
   }
#endif

   copySol(better.Store(), better.Store()+N, sol);
   fitness = fit_better;

   if (eval)
      delete eval;

   delete[] arindex;
   delete[] arindex_raw;

   return counteval;
}

