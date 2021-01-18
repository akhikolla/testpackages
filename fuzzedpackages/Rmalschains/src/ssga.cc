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

#include "ssga.h"
#include <cassert>

using namespace realea;

SSGA::SSGA(Random *random) :  ICrossEAlgorithm(random), m_select(NULL), m_replace(NULL), m_mutation(NULL), 
m_imutation(NULL) {
}

void SSGA::setSelect(ISelect *sel) {
    m_select = sel;
    m_select->setRandom(m_random);
    appendSignal(m_select);
}

void SSGA::setReplacement(IReplace *replace) {
    m_replace = replace;
    appendSignal(m_replace);
}

void SSGA::setMutation(IMutation *mutation) {
    m_imutation = mutation;
 
    if (m_problem) {
    	mutation->setDomain(m_problem->getDomain());
        m_mutation = new Mutation(mutation, m_pmut);
	m_mutation->setRandom(m_random);
	m_mutation->setDomain(m_problem->getDomain());
    }
    
}

void SSGA::setMutationRate(double prob_mutation) {
    m_pmut = prob_mutation;
}

void SSGA::setProblem(Problem *problem) {
    ICrossEAlgorithm::setProblem(problem);
    m_select->setDomain(m_problem->getDomain());

    if (!m_mutation && m_imutation) {
       m_imutation->setDomain(m_problem->getDomain());
       m_mutation = new Mutation(m_imutation);
       m_mutation->setRandom(m_random);
	m_mutation->setDomain(m_problem->getDomain());
    }
}



SSGA::~SSGA(void) {
    if (m_select)
      delete m_select;

    if (m_replace) 
      delete m_replace;

    if (m_mutation) 
      delete m_mutation;
}

unsigned SSGA::init(void) {
    // Primero inicio la población
    m_pop->reset(m_problem->getDomain());
    // Inicio los distintos elementos: Running, Replace, ...
    reset();

   // Compruebo los elementos (mutation no lo compruebo al ser opcional
    if (m_select == NULL) {
	throw new ConfigException("select");
    }
    if (m_replace == NULL) {
	throw new ConfigException("replace");
    }
    if (m_cross == NULL) {
      throw new ConfigException("cross");
    }

    m_pop->eval(m_init_eval);
    return m_running->numEval();
}

unsigned SSGA::realApply(tChromosomeReal &sol, tFitness &fitness) {
    unsigned pos_mom, pos_dad;
    tChromosomeReal crom(m_problem->getDimension());
    tFitness best_fit;
    unsigned better;

    unsigned initMax = m_running->numEval();

   // While is not finish
   while (!m_running->isFinish()) {
        if (m_stat)
	    m_stat->newGeneration();

	// Selecciono los individuos a cruzar
	m_select->select(m_pop, &pos_mom, &pos_dad);
	assert (m_pop->getInd(pos_mom) != m_pop->getInd(pos_dad));

	// Aplico el cruce sobre ambos
	cross(pos_mom, pos_dad, crom);
	// Si es necesario muto
	if (m_mutation) {
	    m_mutation->apply(crom);

	    if (m_stat)
	       m_stat->newEvent("Mutation");
	}

	// Genero el nuevo individuo
	tIndividualRealPtr newind = m_pop->getInstance(crom);
	// Lo evalúo
	m_new_eval->eval(newind);

	// Obtengo el nuevo individuo
	unsigned candreplace = m_replace->getCandidate(m_pop, newind);

	if (m_replace->mustBeReplace(m_pop->getInd(candreplace), newind)) 
	    m_pop->replace(candreplace, newind);
	else 
	    delete newind;

	better = m_pop->getBest();
	best_fit = m_pop->getInd(better)->perf();

	if (m_stat)
	   m_stat->endGeneration(best_fit);

   } // De running

   // Obtengo el mejor
   unsigned pos = m_pop->getBest();
   tIndividualRealPtr best= m_pop->getInd(pos);

   tChromosomeReal bestsol= best->sol();
   copy(bestsol.begin(), bestsol.end(), sol.begin());
   fitness = best->perf();
   unsigned neval = m_running->numEval()-initMax;
   m_running->reset();
   return neval;
}

void SSGA::cross(unsigned pos_mom, unsigned pos_dad, tChromosomeReal &crom) {
    tIndividualReal *mom, *dad;

    mom = m_pop->getInd(pos_mom);
    dad = m_pop->getInd(pos_dad);
    // Llamo al método correspondiente de cruce
    (*m_cross)(mom, dad, crom);
}


