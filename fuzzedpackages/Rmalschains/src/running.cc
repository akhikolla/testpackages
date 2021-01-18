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

#include "running.h"
#include "debug.h"
#include <math.h>
#include <cassert>
#include <cstdio>

using namespace realea;

OptimeCriterion::OptimeCriterion(double optime, double dif) {
      m_optime = optime;
      m_threshold = dif;
}

void OptimeCriterion::setThreshold(double dif) {
      assert(dif >= 0);
      m_threshold = dif;
}

bool OptimeCriterion::isOptime(double fitness) {
	if (fitness < m_optime)
	  return true;
	else 
	    return (fitness-m_optime <= m_threshold);
}

double OptimeCriterion::getThreshold(void) {
	return m_threshold;
}

void Running::increm(void) {
    if (m_optimized) {
      print_error("Warning: Optimized value achieved\n");
//      throw new RunningException("Optimized value achieved");
   }
   if (m_maxmsecs == 0 && m_neval == m_maxeval) {
      print_error("Warning: Max eval achieved\n");
//      throw new RunningException("Max eval achieved");
   }

   m_neval += 1;

}

void Running::reset(void) {
	m_neval = 0;
	m_timeInit = clock();
	m_optimized = false;
}

Running::Running(OptimeCriterion *isOptime) : m_checkOptime(isOptime), m_children() {
	m_maxeval = m_neval = 0;
	m_optimized = false;
	m_parent = NULL;
	m_maxmsecs = 0;
	m_timeInit = 0;
}

void Running::setMaxEval(unsigned int max) {
	m_maxeval = max;
	m_neval = 0;
}

unsigned int Running::maxEval(void) {
	return m_maxeval;
}

bool Running::isOptime(double fitness) {
	if (m_checkOptime->isOptime(fitness)) {
	    m_optimized = true;
	    return true;
	}
	else {
		return false;
	}
	
}

bool Running::isFinish(void) {
	if (m_optimized || (m_maxmsecs == 0 && (m_neval >= m_maxeval)) ) {
	    return true;
	}
	else if (m_parent != NULL) {
	    return m_parent->isFinish();
	}
	else if (m_maxmsecs > 0) {
	    clock_t current = clock();
	    clock_t dif = (10*(current-m_timeInit))/CLOCKS_PER_SEC;
	    return (dif >= m_maxmsecs);
	} 
	else {
	    return false;
	}
}

double Running::ratio(void) {
	if (m_neval == 0) {
		return 0;	
	}
	else {
		return (( (double)m_neval)/m_maxeval);
	}
}

unsigned int Running::numEval(void) {
	return m_neval;
}

void Running::setThreshold(double dif) {

	if (m_neval > 0) {
	    throw new RunningException("Threshold can't be changed in running");
	}
	
	m_checkOptime->setThreshold(dif);	
}

double Running::getThreshold(void) {
	if (m_checkOptime == NULL) {
	    throw new RunningException("Max eval achieved");
	}

	return m_checkOptime->getThreshold();
}

void Running::notifyEval(double fit) {
//    assert(fit >= 0);
    increm();
    if (isOptime(fit))
	m_best = fit;
    else if (m_neval == 1)
      m_best = fit;
    else if (m_checkOptime->isBetter(fit, m_best) ) {
	m_best = fit;
    }

   if (m_parent != NULL) {
      m_parent->notifyEval(fit);
   }

}

RunningPtr Running::getSubRunning(unsigned submaxeval) {
    Running *run = new Running(m_checkOptime);
    run->setMaxEval(submaxeval);
    run->m_parent = this;
    m_children.push_back(run);
    return run;
}


Running::~Running(void) {
   list<Running*>::iterator item;

   for (item = m_children.begin(); item != m_children.end(); ++item) {
	delete *item;
   }
}

void Running::setMaxTime(unsigned ms) {
    m_maxmsecs = ms;
    m_timeInit = clock();
}
