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

#include "localsearch.h"
#include "random.h"
#include <cassert>

using namespace realea::internal;

NewIndividualLocalSearchManager::NewIndividualLocalSearchManager(ILocalSearch *ls) : m_ls(ls) {
}

NewIndividualLocalSearchManager::~NewIndividualLocalSearchManager(void) {
    if (m_ls) {
	delete m_ls;
    }
}

bool RatioLocalSearchManager::applyNewSol(tChromosomeReal &sol, tFitness *pfitness, ILSParameters *par) {
   ILSParameters *params=NULL;
   bool initconf;
   bool mustApply;
   bool applied=false;

    // Check if must be apply the LS
	if (!m_random) {
	    throw new ConfigException("LocalSearch::random");
	}

    mustApply = (m_random->rand() < m_ratio);

    if (mustApply) {
	if (m_ls == NULL) {
	    throw new ConfigException("LocalSearch::ls");
	}

       initconf = (params != NULL);

       if (initconf) {
	  params = m_ls->getInitOptions(sol);
       }
       else {
	  params = par;
       }

       applied = true;
       m_ls->apply(params, sol, *pfitness, m_intensity);

       if (initconf) {
	  delete params;
       }
    }

    return applied;
}

RatioLocalSearchManager::RatioLocalSearchManager(ILocalSearch *ls, unsigned intensity, double ratio) : NewIndividualLocalSearchManager(ls) {
   assert(intensity > 0);
    m_intensity = intensity;
    assert(ratio > 0 && ratio <= 1);
    m_ratio = ratio;
}

unsigned RatioLocalSearchManager::getIntensity(void) {
    return m_intensity;
}

void NewIndividualLocalSearchManager::reset(void) {
}

void RatioLocalSearchManager::setRandom(Random *random) {
	m_random = random;
}

