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

#include "chc.h"

using namespace realea;

CHC::CHC(Random *random): ICrossEAlgorithm(random) {
}

unsigned CHC::getDefaultPopsize(void) {
	return 50;
}


CHC::~CHC(void) {
}


unsigned CHC::init(void) {
  // Inicio la población
   m_pop->reset(m_problem->getDomain());

    // Inicio los distintos elementos: Running, Replace, ...
    reset();

   if (m_cross == NULL) {
      throw new ConfigException("cross");
   }

   m_running->reset();
   m_pop->eval(m_init_eval);

   m_init_threshold = m_problem->getDimension()*BITS_GEN/4;
   m_threshold = m_init_threshold;
   return m_running->numEval();
}


unsigned int CHC::realApply(tChromosomeReal &sol, tFitness &fitness) {
   tFitness fitbest;

   m_running->reset();

  // While is not finish
   while (!m_running->isFinish()) {
	m_pop->random();

	if (m_stat) 
	   m_stat->newGeneration();

	// Cruzo (añado los nuevos al final)
	int crossed = cross(m_threshold);
	// Evalúo los últimos
	m_pop->eval(m_new_eval, m_running->pending());

	// Elimino los peores
	m_pop->removeWorses();

	fitbest = m_pop->getInd(0)->perf();

	if (m_stat)
	   m_stat->endGeneration(fitbest);

	// Compruebo si debe de reiniciar
	if (crossed == 0 && !m_running->isFinish()) {
	    m_threshold -= 1;

	    if (m_threshold < 0) {
		((PopulationRealCHC *)m_pop)->restart(m_problem->getDomain());
		m_threshold = m_init_threshold; 
		m_pop->eval(m_init_eval);
		if (m_stat) 
		   m_stat->newEvent("Restart");
	    }
	} // De comprobar si debe de cruzar

   } // De running

   // Obtengo el mejor
   unsigned pos = m_pop->getBest();
   tIndividualRealPtr best= m_pop->getInd(pos);

   tChromosomeReal bestsol= best->sol();
   copy(bestsol.begin(), bestsol.end(), sol.begin());
   fitness = best->perf();
   return m_running->numEval();
} 

unsigned int CHC::distHamming(tIndividualRealPtr ind1, tIndividualRealPtr ind2) {
    tIndividualRealCHC* mom, *dad;
    DomainRealPtr domain = m_problem->getDomain();

    mom = (tIndividualRealCHC*)ind1;
    dad = (tIndividualRealCHC*)ind2;
    mom->calculateBin(domain);
    dad->calculateBin(domain);
 
    return mom->distHammingOpt(dad);
}

unsigned CHC::cross(unsigned threshold) {
    int dad, mom;
    unsigned size = m_pop->size();
    unsigned count=0;

    for (unsigned i = 0; i < size/2; i++) {
	mom = 2*i;
	dad = 2*i+1;

	unsigned dist = distHamming(m_pop->getInd(mom), m_pop->getInd(dad));

	if (dist > threshold*2) {
	    tChromosomeReal sol(m_problem->getDimension());
	    (*m_cross)(m_pop->getInd(mom), m_pop->getInd(dad), sol);
	    tIndividualReal *ind = m_pop->getInstance(sol);
	    m_pop->append(ind);
	    ind->setId(m_pop->size());
	    count++;
	}
    }

    return count;
}

/**
 * Para CHC esta función no tiene sentido, ya que será siempre de 50 individuos
 */
void CHC::setPopsize(unsigned int popsize) {
    if (m_pop) {
	delete m_pop;
    }

    m_pop = new PopulationRealCHC(m_random,2*popsize, popsize);
}
