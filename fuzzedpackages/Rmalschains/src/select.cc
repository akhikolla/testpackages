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

#include "select.h"
#include "random.h"
#include "distance.h"
#include <cassert>

using namespace realea;

SelectNAM::SelectNAM(unsigned nam) : m_num(nam) {
    assert(m_num > 0);
}

void SelectNAM::select(PopulationRealPtr pop, unsigned *pmom, unsigned *pdad) {
    int max = pop->size();
    int ndim = pop->ndim();
    tIndividualReal *mom, *elem;
    unsigned pos_mom, pos_elem, pos_best=0;
    double dist, best_dist=0;
    bool *checkGen = new bool[ndim];
    
    m_domain->getSearchDomain(checkGen, ndim);
    fill(checkGen, checkGen+ndim, true);

    assert(m_num > 0);
    int *sample = new int[max];
    // Obtengo un elemento aleatoriamente
    initSample(sample, max);
    pos_mom = m_random->getSample(sample, &max);
    mom = pop->getInd(pos_mom);

    for (unsigned i = 0; i < m_num; ++i) {
	pos_elem = m_random->getSample(sample, &max);
	elem = pop->getInd(pos_elem);
	dist = distreal(mom->sol(), elem->sol(), checkGen);

	if (i == 0 || (dist > best_dist) ) {
	    pos_best = pos_elem;
	    best_dist = dist;
	}
    } // De recorrer la población

    *pmom = pos_mom;
    *pdad = pos_best;
    delete[] checkGen;
    delete[] sample;
}

SelectTournament::SelectTournament(unsigned nam) : m_num(nam) {
    assert(m_num > 0);
}

/**
 * Apply a tournament, returning the individual with better fitness
 *
 * @param pop population
 * @param tournament_size tournament size
 * @param random random generator
 * @param positions positions of individuals
 * @param maxpos position's size 
 *
 * @return the best individual
 */
tIndividualRealPtr applyTournament(PopulationRealPtr pop, unsigned tournament_size, Random *random, int *positions, int *maxpos) {
    tIndividualRealPtr ind=NULL;

    int posi = random->getSample(positions, maxpos);
    tIndividualRealPtr best = pop->getInd(posi);

    for (unsigned i = 1; i < tournament_size; i++) {
	posi = random->getSample(positions, maxpos);
	ind = pop->getInd(posi);

	if (ind->isBetter(best)) {
	   best = ind;
	}
    }

    return ind;
}

void SelectTournament::select(PopulationRealPtr pop, unsigned *pmom, unsigned *pdad) {
    int max = pop->size();
    int *positions = new int[max];
    tIndividualReal *mom, *dad;
    initSample(positions, max);

    mom = applyTournament(pop, m_num, m_random, positions, &max);
    dad = applyTournament(pop, m_num, m_random, positions, &max);

    *pmom = mom->getId();
    *pdad = dad->getId();
    delete[] positions;
}
