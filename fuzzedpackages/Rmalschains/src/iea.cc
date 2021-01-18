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

#include "iea.h"
#include <cassert>

using namespace realea;
using namespace realea::internal;

ClassEAlgorithm::ClassEAlgorithm(Random *random) : m_problem(NULL), m_init_eval(NULL), m_new_eval(NULL), m_running(NULL), m_pop(NULL), m_stat(NULL) {
	setRandom(random);
	m_maxeval = 0;
}

void ClassEAlgorithm::setNewEval(IEval *eval) {
   assert(m_running != NULL);
   assert(m_new_eval == NULL);
   m_new_eval = new EvalRunning(eval, m_running);    
}

void ClassEAlgorithm::setInitEval(IEval *eval) {
   assert(m_running != NULL);
   assert(m_init_eval == NULL);
   m_init_eval = new EvalRunning(eval, m_running);    
}

void ClassEAlgorithm::setMaxEval(unsigned int maxeval) {
	m_maxeval = maxeval;

	if (m_running) {
	   m_running->setMaxEval(m_maxeval);
	}
}

unsigned int ClassEAlgorithm::getMaxEval(void) {
	return m_maxeval;
}


ClassEAlgorithm::~ClassEAlgorithm(void) {
    if (m_pop) 
      delete m_pop;

    if (m_init_eval) 
      delete m_init_eval;

    if (m_new_eval)
      delete m_new_eval;

    m_new_eval = m_init_eval = NULL;
}

void ClassEAlgorithm::setPopsize(unsigned int popsize) {
    if (m_pop) {
	delete m_pop;
    }

    m_pop = new PopulationReal(m_random, popsize, popsize);
}

void ClassEAlgorithm::reset(void) {
    if (m_running)
       m_running->reset();
}

void ClassEAlgorithm::setRunning(Running *running) {
	if (m_running != NULL) {
	    delete m_running;
	}
    
	m_running = running;

	if (m_maxeval > 0 && m_running != NULL) {
	   m_running->setMaxEval(m_maxeval);
	}

}

RunningPtr ClassEAlgorithm::getRunning(void) {
     assert(m_running != NULL);
     return m_running;
}

void ClassEAlgorithm::storeIndividual(tIndividualRealPtr ind, tGen **paind, unsigned *pmax, 
					tGen **paoptional, unsigned *psize_optional) {
    vector<tGen> sol = ind->sol();
    unsigned size = sol.size()+1;

    tGen *aind = new tGen[size];
    copy(sol.begin(), sol.end(), aind);
    aind[size-1] = ind->perf();
    
    *pmax = size;
    *paind = aind;
    *paoptional = NULL;
    *psize_optional= 0;
}

void ClassEAlgorithm::recoverIndividual(unsigned pos, tGen *aind, unsigned size, tGen *aoptional, unsigned optional_size) {
    assert(size > 1 && (size == m_pop->ndim()+1) );
    vector<tGen> sol(size-1);
 
    copy(aind, aind+size-2, sol.begin());
    tFitness fitness = aind[size-1];

    tIndividualRealPtr ind = m_pop->getInstance(sol, fitness);
    m_pop->replace(pos, ind);
}

ICrossEAlgorithm::ICrossEAlgorithm(Random *random) : ClassEAlgorithm(random), m_cross(NULL), m_icross(NULL) {
}

ICrossEAlgorithm::~ICrossEAlgorithm(void) {
   if (m_cross) 
      delete m_cross;
}

void ICrossEAlgorithm::setCross(ICrossBinaryPtr cross) {
    m_icross = cross;
}

void ICrossEAlgorithm::setProblem(Problem *problem) {
    ClassEAlgorithm::setProblem(problem);
}


void ICrossEAlgorithm::reset(void) {
    if (m_cross) {
       m_cross->reset();
    }
    else if (m_icross) {
       m_icross->setRandom(m_random);	
       m_icross->setDomain(m_problem->getDomain());
       m_icross->setRunning(m_running);
       m_cross = new CrossBinary(m_icross);
       appendSignal(m_cross);
    }
    
    ClassEAlgorithm::reset();	    
}
