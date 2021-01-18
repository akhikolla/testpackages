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
#include "simplex.h"
#include "distance.h"
#include "random.h"
#include <functional>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <cmath>

using namespace realea;

class SimplexParams : public ILSParameters {
public:
    /**
     * Store the solution over whose the different operators will be applied
     */
    vector<tChromosomeReal> m_simplex;

    /**
     * Store for each solution in m_simplex its fitness
     */
    vector<tFitness > m_fitnessElem;

    /**
     * Store for each dimension the accumulated distances stores in m_simplex
     */
    vector<tGen> m_cum_simplex;

    int ndim(void) {
        assert(!m_simplex.empty());
	return m_simplex[0].size();
    }
    
    SimplexParams(void) : m_simplex(), m_fitnessElem(), m_cum_simplex() {
    }

    void calculateCentroide(void);
    void getBest(vector<tGen> &sol, tFitness &fitness); 
};

/**
 * Let sort a vector index in function of an indicated vector (passed by the constructor)
 */
template<class T>
class SortIndex : binary_function<int,int,bool>{
public:
    T* m_values;

    SortIndex(T* values) : m_values(values) {}

    bool operator()(int &elem1, int &elem2) {
	return m_values[elem1] < m_values[elem2];
    }
};


/**
 * Esta función permite realiza la suma acumulada de cada solución
 *
 * @param result resultado acumulado
 * @param elem nuevo vector solución a sumar
 */
tChromosomeReal &accumulate_var(tChromosomeReal &result, tChromosomeReal &elem) {
   int i, ndim=result.size();

   for (i = 0; i < ndim; ++i) {
      result[i] += elem[i];
   }

   return result;
}



void SimplexParams::calculateCentroide(void) {
   int ndim = m_simplex[0].size();

   if (m_cum_simplex.empty()) {
      fill_n(back_inserter(m_cum_simplex), ndim, 0.0);
   }
   else {
      // Inicio el centroide
   fill(m_cum_simplex.begin(), m_cum_simplex.end(), 0.0);
   }
     
   // acumulo los valores por cada variable
   m_cum_simplex = accumulate(m_simplex.begin(), m_simplex.end(), m_cum_simplex, accumulate_var);
}

void SimplexParams::getBest(vector<tGen> &sol, tFitness &fitness) {
   int i, ndim=sol.size();
   int best;
   // Obtengo los índices
   vector<int> index(ndim+1);

   for (i = 0; i <= ndim; i++) {
      index[i] = i;
   }

   SortIndex<tFitness > mySort(&m_fitnessElem[0]);
   
   // obtengo el mejor vértice
   best = *min_element(index.begin(), index.end(), mySort);

   // Realizo la copia
   copy(m_simplex[best].begin(), m_simplex[best].end(), sol.begin());
   fitness = m_fitnessElem[best];
}

Simplex::Simplex(void) {
}

void Simplex::getExtremes(ILSParameters *params, int &best_index, int &next_worst_index, int &worst_index) {
   SimplexParams *p = (SimplexParams *) params;
   tFitness bestFitness, worstFitness, nextworstFitness;
   tFitness fitness;
   int elem = p->m_fitnessElem.size();
   int i;

   best_index = worst_index = next_worst_index = 0;
   bestFitness = worstFitness = nextworstFitness = p->m_fitnessElem[0];

   for (i = 1; i < elem; i++) {
      fitness = p->m_fitnessElem[i];

      if (m_problem->isBetter(fitness, bestFitness)) {
	 best_index = i;
	 bestFitness = fitness;
      }
      else if (m_problem->isBetter(worstFitness, fitness)) {
	 next_worst_index = worst_index;
	 nextworstFitness = worstFitness;
	 worst_index = i;
	 worstFitness = fitness;
      }
      else if (m_problem->isBetter(nextworstFitness, fitness) && i != worst_index) {
	 next_worst_index = i;
	 nextworstFitness = fitness;
      }

   }

}

ILSParameters *Simplex::getInitOptions(tChromosomeReal &sol) {
    // Copio la solución actual
   int ndim;

   ndim = sol.size();

   SimplexParams *p = new SimplexParams();

   p->m_simplex.clear();
   p->m_simplex.reserve(ndim+1);
   p->m_fitnessElem.reserve(ndim+1);

   return p;
}

tFitness Simplex::move(ILSParameters *params, int posi, double factor) {
   SimplexParams *p = (SimplexParams *) params;
   tGen *worst_point;
   tFitness *worst_value;
   unsigned ndim = p->ndim();
   tFitness rfac,new_value;
   /* Generate a new point */
   vector<tGen> new_point(ndim);
   unsigned i;

   worst_point = &p->m_simplex[posi][0];
   worst_value = &p->m_fitnessElem[posi];

   rfac = (1.0-factor)/(double)ndim;
   DomainRealPtr domain = m_problem->getDomain();

   for (i=0;i<ndim; i++) {
      // Utilizo factor y el centroide para moverme
      new_point[i]=p->m_cum_simplex[i]*rfac - worst_point[i]*(rfac-factor);
      new_point[i]=domain->clip(i, new_point[i]);
   }

   new_value = m_eval->eval(new_point);

   /* Check worst_point replacement */
   if (m_problem->isBetter(new_value, *worst_value)) {
      *worst_value = new_value;

      // Actualizo el centroide y el peor elemento
      for (i=0; i<ndim; i++) {
	 p->m_cum_simplex[i] += new_point[i]-worst_point[i];
	 worst_point[i] = new_point[i];
      }
   }

   return new_value;
}

unsigned Simplex::restart_simplex(ILSParameters *params, int best, unsigned max) {
   SimplexParams *p = (SimplexParams *) params;
   int ndim = p->ndim();
   int elem = p->m_simplex.size();
   int i, j;
   unsigned num=0;

   for (i = 0; i < elem && num < max && !m_running->isFinish(); ++i) {
      if (i == best)
	continue;

      // Lo acerca al mejor
      for (j = 0; j < ndim; ++j) {
	 p->m_simplex[i][j] = 0.5*(p->m_simplex[i][j]+p->m_simplex[best][j]);
      }

      p->m_fitnessElem[i] = m_eval->eval(p->m_simplex[i]);
      num++;
   }
   
   return num;
}


unsigned Simplex::apply(ILSParameters *params, tChromosomeReal &sol, tFitness &sol_perf, unsigned maxeval) {
   SimplexParams *p = (SimplexParams*) params;
   tFitness newFitness;
   int best, worst, next_worst;
   double factor = 1.0;
   unsigned num;

   num = 0;

   if (p->m_simplex.empty()) {
      num = initParams(sol, sol_perf, p);
   }

   p->calculateCentroide();

   while (num < maxeval && !m_running->isFinish()) {
       getExtremes(p, best, next_worst, worst);

       tFitness &bestFitness = p->m_fitnessElem[best];
       tFitness &worstFitness = p->m_fitnessElem[worst];
       tFitness &nextworstFitness = p->m_fitnessElem[next_worst];

       int actual = worst;

       // Modifico la peor solución actual mediante reflexión
       newFitness = move(p, actual, -factor);
       num++;

       // Si mejora al mejor aplica expansión
       if (m_running->isFinish()) {

       }
       else if (m_problem->isBetter(newFitness, bestFitness) || newFitness == bestFitness) {
	  newFitness = move(p, actual, 2.0*factor);
	  num++;
       }
       // Si es peor que el segundo peor de la población modifica la escala
       else if (!m_problem->isBetter(newFitness, nextworstFitness)) {
	  // Hace contracción
	  newFitness = move(p, actual, factor/2.0);
	  num++;

	  if (m_problem->isBetter(worstFitness, newFitness) || worstFitness == newFitness) {
	     // Reinicia simplex
	     num += restart_simplex(p, best, maxeval-num);
	     // Recalculo el centroide
	     p->calculateCentroide();
	  }
       }
   }

    p->getBest(sol, sol_perf);
    return num;
}

unsigned SimplexDim::initParams(tChromosomeReal &sol, tFitness fitness, ILSParameters *params) {
   SimplexParams *p = (SimplexParams *) params;
   vector<tGen> copysol(sol);
   tFitness newFitness;
   tGen min, max, dif;
   int ndim, i;

   // Almaceno la solución
   p->m_simplex.push_back(sol);
   p->m_fitnessElem.push_back(fitness);

   DomainRealPtr domain = m_problem->getDomain();
   ndim = domain->getDimension();
     
   // Guarda ndim variaciones del vector como conjunto 
   for (i = 0; i < ndim; i++) {
      domain->getValues(i, &min, &max);
      dif = (max-min);
      copysol[i] += 0.1*dif;
      copysol[i] = domain->clip(i, copysol[i]);
      newFitness = m_eval->eval(copysol);
      p->m_fitnessElem.push_back(newFitness);
      p->m_simplex.push_back(copysol);
      copysol[i] = sol[i];
   }

   return ndim;
}

unsigned SimplexNeigh::initParams(tChromosomeReal &sol, tFitness fitness, ILSParameters *params) {
   SimplexParams *p = (SimplexParams *) params;
   vector<tGen> copysol(sol);
   tFitness newFitness;
   int ndim, i;

   // Almaceno la solución
   p->m_simplex.push_back(sol);
   p->m_fitnessElem.push_back(fitness);

   DomainRealPtr domain = m_problem->getDomain();
   ndim = domain->getDimension();

   // Guardo ndim soluciones más cercanas por dimensión
   vector<unsigned> minind(sol.size());

   min_dim_distance(sol, m_pop, minind);

   for (i = 0; i < ndim; i++) {
      unsigned pos = minind[i];
      tChromosomeReal ind = m_pop->getInd(pos)->sol();
      copy(ind.begin(), ind.end(), copysol.begin());
      newFitness = m_eval->eval(copysol);
      p->m_fitnessElem.push_back(newFitness);
      p->m_simplex.push_back(copysol);
   }

   return ndim;

}


