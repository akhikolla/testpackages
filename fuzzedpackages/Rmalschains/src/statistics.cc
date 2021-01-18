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

#include "statistics.h"
#include "debug.h"
#include <cassert>
#include <cstdio>

using namespace realea;

Statistics::Statistics(int num_gen) {
    reset();
    m_rate = num_gen;
}

void Statistics::reset(void) {
    m_experiment = m_generation = 0;
}

void Statistics::newExperiment(void) {
    print_info("Experiment: %d\n", m_experiment);
    ++m_experiment;
    m_generation = 0;
}

void Statistics::newGeneration(void) {
    ++m_generation;
}

void Statistics::endGeneration(tFitness best) {
    if (m_generation > 1) {
        if (m_problem->isBetter(m_lastbest, best)) {
            print_info("m_lastBest: %Le\tbest : %Le\n", m_lastbest, best);
	}

	assert(!m_problem->isBetter(m_lastbest, best));
    }

    if (m_rate > 0) {
	unsigned rate = (unsigned) m_rate;

	if ((m_generation % rate)==0) {
	   print_info("Best[%d]: %Le\n", m_generation, best);
	}

    }

    m_lastbest = best;
}

void Statistics::newEvent(string event) {
    map<string,bool>::iterator pos;
    bool active;

    pos = m_events.find(event);

    if (pos == m_events.end()) 
      active = false;
    else {
	active = (*pos).second;	
    }

    if (active) 
	print_info("%s:[%d]\n", event.c_str(), m_generation);
}

void Statistics::endExperiment(void) {
    print_info("BestExperiment[%d]: %Le\n", m_experiment, m_lastbest);
}

void Statistics::activeEvent(string event) {
    m_events[event] = true;
}
