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

#include "problem.h"
#include <cassert>
#include "running.h"

using namespace realea;

Problem::Problem(void) : m_checkOptime(NULL), m_domain(NULL), m_notify(NULL) {
	
}

bool Problem::isBetter(tFitness value1, tFitness value2) {
	return m_checkOptime->isBetter(value1, value2);
}

bool Problem::minimize(void) {
   if (m_domain == NULL) {
      throw new ConfigException("domain");
   }
 
    return m_checkOptime->minimize();
}

bool Problem::maximize(void) {
    return !minimize();
}

void Problem::setMaximize(void) {
    m_checkOptime->setMaximize();
}

void Problem::setMinimize(void) {
    m_checkOptime->setMinimize();
}


void Problem::setDimension(unsigned int dim) {
    m_domain = new DomainReal(dim);
}

void Problem::setDomainValues(unsigned int gen, tGen min, tGen max, bool check) {
   if (m_domain == NULL) {
      throw new ConfigException("domain");
   }

   m_domain->setValues(gen, min, max, check);
}


Problem::~Problem(void) {
	if (m_checkOptime) {
		delete m_checkOptime;
	}

	if (m_domain) {
	    delete m_domain;
	}

}

void Problem::setOptimize(tFitness optime, double threshold) {
	m_checkOptime = new OptimeCriterion(optime,threshold);
}

DomainRealPtr Problem::getDomain(void) {
    if (m_domain == NULL) {
	throw new ConfigException("domain");
    }
    return m_domain;
}

void Problem::setThreshold(double dif) {
    if (m_checkOptime == NULL) {
	throw new ConfigException("optime");
    }

    m_checkOptime->setThreshold(dif);
}

tFitness Problem::getOptime(void) {
    if (m_checkOptime == NULL) {
	throw new ConfigException("optime");
    }
    return m_checkOptime->getOptime();
}

void Problem::setEval(tEval eval) {
    m_eval = eval;
}

void Problem::setMaxEval(unsigned eval) {
    m_maxeval = eval;
}

unsigned Problem::getMaxEval(void) {
	return m_maxeval;
}


void Problem::setAfterEval(tAfterEvalFunction function) {
	m_notify = function;
}


tFitness Problem::eval(const double *sol, int dim) {
    tFitness fit = (*m_eval)(sol, dim);

    if (m_notify) {
	(*m_notify)(sol, dim, fit);
    }
    return fit;
}


tFitness Problem::eval(const tChromosomeReal &sol) {
    tFitness fit = (*m_eval)(&sol[0], sol.size());

    if (m_notify) {
	(*m_notify)(&sol[0], sol.size(), fit);
    }
 
    return fit;
}

unsigned Problem::getDimension(void) {
    return m_domain->getDimension();
}

void Problem::copy(Problem *problem) {
   assert(problem->m_checkOptime ==  NULL &&
	  problem->m_domain == NULL);
    
   problem->m_checkOptime = this->m_checkOptime;
   problem->m_domain = this->m_domain;
   problem->m_eval = this->m_eval;
   problem->m_maxeval = this->m_maxeval;
}
