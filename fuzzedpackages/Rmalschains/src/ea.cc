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

#include "ea.h"
#include <cassert>

using namespace realea;
using namespace realea::internal;




EAlgorithm::EAlgorithm(IEAlgorithm *alg, ProblemParamPtr problem) : m_alg(alg), default_popsize(alg->getDefaultPopsize()), m_cross(NULL), m_stat(NULL) {
	setProblem(problem);
}

void EAlgorithm::setShow(Statistics *stat) {
    if (m_problem) 
       stat->setProblem(m_problem);

    m_alg->setStat(stat);  
    m_stat = stat;
}

void EAlgorithm::setEval(IEval *eval) {
    m_eval = eval;
}

void EAlgorithm::setProblem(ProblemParamPtr problem) {
     m_problem = problem.get();
     m_alg->setProblem(m_problem);
     m_running = new Running(m_problem->getFinishCriterion());
     m_running->setMaxEval(m_problem->getMaxEval());
     m_alg->setRunning(m_running);
     m_alg->setMaxEval(m_running->maxEval());
     m_alg->setInitEval(m_problem);
     m_alg->setNewEval(m_problem);

    // Set the criterion to the individual class
    if (problem->minimize())
	tIndividualReal::setMinimize();
    else 
	tIndividualReal::setMinimize();

    appendSignal(m_alg);

    if (m_stat) 
       m_stat->setProblem(m_problem);
}

EAlgorithm::~EAlgorithm(void) {

    if (m_alg)
       delete m_alg;

    if (m_stat)
      delete m_stat;

    if (m_running) 
      delete m_running;
}

void EAlgorithm::setMaxEval(unsigned int maxeval) {
   m_alg->setMaxEval(maxeval);
}

void EAlgorithm::setPopsize(unsigned int popsize) {
    m_alg->setPopsize(popsize);
}

RunningPtr EAlgorithm::getRunning(void) {
    return m_running;
}

void EAlgorithm::setDefaultPopsize(void) {
   unsigned maxeval = m_alg->getMaxEval();
   unsigned max;

   if (maxeval < default_popsize) {
      max = maxeval;
   }
   else {
      max = default_popsize;
   }

   // Asigna la población con el tamaño indicado
   m_alg->setPopsize(max);
}

/**
 * Este método es el llamado para resolver el algoritmo evolutivo, si no se especificó una llamada de popsize
 * le asigna el menor valor
 *
 */
unsigned EAlgorithm::apply(tChromosomeReal &sol, tFitness *pfitness) {
   // Si no se ha definido la población lo define
   if (m_alg->getPop() == NULL) 
      setDefaultPopsize();

   if (m_problem == NULL) {
      throw new ConfigException("problem");
   }

   if (m_stat)
       m_stat->newExperiment();

   unsigned res = m_alg->init();
   
   res += m_alg->realApply(sol, *pfitness);

   if (m_stat)
     m_stat->endExperiment();

   return res;
}
