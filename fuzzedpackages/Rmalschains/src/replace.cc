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

#include "replace.h"
#include "distance.h"
#include "debug.h"
#include <cassert>
#include <cstdio>

using namespace realea;

IReplace::IReplace(void) : m_total(0), m_success(0) {

}

void IReplace::reset(void) {
    double ratio;

    if (m_total != 0) {
      ratio = ((double)m_success)/m_total;
      print_info("Replacement Ratio: %.0lf%%\n", 100*ratio);
    }

    m_total = m_success = 0;
    print_info("IReplace: Se reinicia\n");
}

bool IReplace::mustBeReplace(tIndividualRealPtr old, tIndividualRealPtr newind) {
    bool changed = newind->isBetter(old);

    if (changed) {
	m_success += 1;
    }

    m_total += 1;
    return changed;
}

unsigned ReplaceWorst::getCandidate(PopulationRealPtr pop, tIndividualRealPtr newind) {
    return pop->getWorst();
}

unsigned ReplaceDC::getCandidate(PopulationRealPtr pop, tIndividualRealPtr newind) {
    unsigned posind;
    distanceMin(newind->sol(), pop, &posind);
    assert(posind < pop->size());
    return posind;
}
