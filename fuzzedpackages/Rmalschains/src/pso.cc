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
#include "pso.h"
#include "populationpso.h"
#include <cassert>

using namespace realea;

void PSO::setConfigPSO(ConfigPSO *config) {
    setGlobalConfigPSO(config);
    m_config = config;
}

PSO::~PSO(void) {
   delGlobalConfigPSO();
}

PSO::PSO(Random *random) : ClassEAlgorithm(random),m_config(NULL) {
}

unsigned PSO::getDefaultPopsize(void) {
    return 30;
}

unsigned PSO::init(void) {
    if (m_config == NULL) {
	setConfigPSO(new ConfigPSO(m_problem->getDomain(), 0.4, 0.98));
    }
    
    // Init the population
    m_pop->reset(m_problem->getDomain());
    m_running->reset();
    m_pop->eval(m_init_eval);
    return m_running->numEval();
}

void PSO::setPopsize(unsigned int popsize) {
    if (m_pop) {
	delete m_pop;
	m_pop = NULL;
    }

    m_pop = new PopulationPSO(m_random, popsize, popsize);
}

unsigned PSO::realApply(tChromosomeReal &sol, tFitness &fitness) {
    // Move all the individuals
    PopulationPSO *pop = (PopulationPSO *) m_pop;

    m_running->reset();

   // While is not finish
   while (!m_running->isFinish()) {
      // Move all the population
      pop->move(m_new_eval, m_running);
   }

   // Obtengo el mejor
   unsigned pos = m_pop->getBest();
   tIndividualRealPtr indbest= m_pop->getInd(pos);

   tChromosomeReal bestsol= indbest->sol();
   copy(bestsol.begin(), bestsol.end(), sol.begin());
   fitness = indbest->perf();

   return m_running->numEval();
}

