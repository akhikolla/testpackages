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
#include "reducepopulation.h"
#include <cassert>

using namespace realea;

PopulationReductionStrategy::PopulationReductionStrategy(void) {
    m_running = NULL;
    m_maxReductions = 0;
    m_numReductions = 1;
}

void PopulationReductionStrategy::setNumReductions(int numReductions) {
    m_maxReductions = numReductions;
}

void PopulationReductionStrategy::config(Running *running) {
    m_running = running;
    assert(m_running != NULL);

    m_maxFES = m_running->maxEval();
    m_initEvals = m_running->numEval();
    m_countdown = m_maxFES/(m_maxReductions+1);
}

bool PopulationReductionStrategy::updatePopulationSize(PopulationReal *pop) {
    int NP = pop->size();
    int evals = m_running->numEval();

    /*if :
     *  -   m_running exists
     *  -   the reductions are to be performed
     *  -   the size of the pop is bigger than 10
     *  -   the number of evaluation has reached the threshhold (m_countdown)
     *  -   all the reductions are not perormed yet
     */

    if (m_running == NULL || m_maxReductions == 0 || NP <= 10
            || evals<=m_countdown ||m_numReductions > m_maxReductions) {
	return false;
    }

    //increment the number of reduction
    m_numReductions++;
    //calculate the new threshold to reach before next reduction
    m_countdown = m_numReductions*m_maxFES/(m_maxReductions+1);

   // *** new reduction ***
   pop->reduceHalf();
   
   return true;
}
