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
#include "selectlsimp.h"
#include "distance.h"
#include <algorithm>
#include <cassert>

using namespace realea;

bool isBetter(tIndividualReal *a, tIndividualReal *b) {
   if (a->isBetter(b)) {
      return true;
   }
   else {
      return false;
   }
}

unsigned SelectBestToImprove::selectIndToImprove(deque<tIndividualReal*> &individuals) {
    deque<tIndividualReal*>::iterator elem;

    assert(individuals.size()>0);

    elem = min_element(individuals.begin(), individuals.end(), isBetter);

    return ( (*elem)->getId());
}

void SelectBestToImprove::getIndsToImprove(PopulationReal *pop, deque<tIndividualReal*> &subpop) {
    unsigned popsize = pop->size();
    tIndividualReal* ind;
    unsigned i, max;

    subpop.clear();

    // It select all individuals that have not been applied the LS and not improve
    for(i = 0, max=0; i < popsize; i++) {
	ind = pop->getInd(i);

	if (ind->getCount("non_improved")==0 && ind->isEval()) {
	    subpop.push_back(ind);
	    max++;
	}
    } 
}

void SelectWithDiversityToImprove::getIndsToImprove(PopulationReal *pop, deque<tIndividualReal*> &subpop) {
    unsigned popsize = pop->size();
    tIndividualReal* ind;
    unsigned i, max;

    subpop.clear();
    m_previous.clear();

    // It select all individuals that have not been applied the LS and not improve
    for(i = 0, max=0; i < popsize; i++) {
	ind = pop->getInd(i);

	if (!ind->isEval()) {
	    continue;
	}

	if (ind->getCount("non_improved")==0) {
	    subpop.push_back(ind);
	    max++;
	}
	else {
	    m_previous.push_back(ind);
	}
    }
}


class GetDistant {
   public:
	GetDistant(deque<tIndividualReal*> references) : m_references(references) {}
	bool operator()(tIndividualReal *a, tIndividualReal *b) {
	    double distA, distB;

	    distA = minDistance(a);
	    distB = minDistance(b);

	    if (distA > distB) {
		return true;
	    }
	    else if (distA < distB) {
		return false;
	    }
	    else {
		return isBetter(a, b);
	    }
	}

	double minDistance(tIndividualReal *ind);

    private:
	deque<tIndividualReal *> m_references;
};

double GetDistant::minDistance(tIndividualReal *ind) {
    double dist, distMin=-1;
    deque<tIndividualReal*>::iterator pos_ref;
    tIndividualReal *ind_ref;

    for (pos_ref = m_references.begin(); pos_ref != m_references.end(); ++pos_ref) {
        ind_ref = (*pos_ref);

	if (ind->isEval()==false) {
	    continue;
	}

        if (ind_ref->perf() == ind->perf()) 
	   continue;

	dist = distreal(ind_ref->sol(), ind->sol());

	if (distMin < 0 || dist < distMin) {
	    distMin = dist;
	}

    }

    return distMin;
}

bool isImproving(tIndividualReal *ind) {
    return (ind->getCount("ls") > 0);
}

deque<tIndividualReal*>::iterator more_distant(deque<tIndividualReal*> &individuals, GetDistant &distant) {
    deque<tIndividualReal*>::iterator posi, greater;
    double dist, distMax;

    if (individuals.empty()) 
       return individuals.end();

    posi = individuals.begin();
    greater = posi;
    distMax = distant.minDistance(*greater);
    posi++;

    while (posi != individuals.end()) {
	dist = distant.minDistance(*posi);

	if (dist > distMax) {
	    greater = posi;
	    distMax = dist;
	}
	posi++;
    }

    return greater;
}

unsigned SelectDiverseToImprove::selectIndToImprove(deque<tIndividualReal*> &individuals) {
    deque<tIndividualReal*>::iterator elem;

    assert(individuals.size()>0);

    if (m_previous.empty()) {
	elem = min_element(individuals.begin(), individuals.end(), isBetter);
    }
    else {
        elem = find_if(individuals.begin(), individuals.end(), isImproving);
    
	if (elem == individuals.end()) {
	    GetDistant distant(m_previous);
	    elem = more_distant(individuals, distant);
	}
    }

    return ( (*elem)->getId());
}

class SortInd{
    public:
	bool operator()(tIndividualReal *a, tIndividualReal *b) {
	    if (a->isEval() && b->isEval())
		return (a->isBetter(b));
	    else if (a->isEval()) 
		return true;
	    else
		return false;
	}
};


unsigned SelectDistantBestToImprove::selectIndToImprove(deque<tIndividualReal*> &individuals) {
    deque<tIndividualReal*>::iterator elem;
    unsigned elem_id;

    assert(individuals.size()>0);

    if (m_previous.empty()) {
	elem = min_element(individuals.begin(), individuals.end(), isBetter);
	elem_id = (*elem)->getId();
    }
    else {
        elem = find_if(individuals.begin(), individuals.end(), isImproving);
	deque<tIndividualReal*>::iterator posi_max;

	if (elem != individuals.end()) {
	    elem_id = (*elem)->getId();
	}
	else {
	   vector<tIndividualReal*> list_bests(individuals.size());
	   vector<tIndividualReal*>::iterator elem_dist;
	   unsigned max = m_maxbests;

	   if (max > individuals.size()) {
	      max = individuals.size();
	   }

	   copy(individuals.begin(), individuals.end(), list_bests.begin());
	   partial_sort(list_bests.begin(), list_bests.begin()+max, list_bests.end(), SortInd());
	   GetDistant distant(m_previous);
	   elem_dist = min_element(list_bests.begin(), list_bests.begin()+max, distant);
	   elem_id = (*elem_dist)->getId();
	}
    }

    return elem_id;
}

unsigned SelectBestDistantToImprove::selectIndToImprove(deque<tIndividualReal*> &individuals) {
    deque<tIndividualReal*>::iterator elem;
    unsigned elem_id;

    assert(individuals.size()>0);

    if (m_previous.empty()) {
	elem = min_element(individuals.begin(), individuals.end(), isBetter);
	elem_id = (*elem)->getId();
    }
    else {
        elem = find_if(individuals.begin(), individuals.end(), isImproving);
	deque<tIndividualReal*>::iterator posi_max;

	if (elem != individuals.end()) {
	    elem_id = (*elem)->getId();
	}
	else {
	   vector<tIndividualReal*> list_bests(individuals.size());
	   vector<tIndividualReal*>::iterator elem_dist;
	   unsigned max = m_maxdists;

	   if (max > individuals.size()) {
	      max = individuals.size();
	   }


	   copy(individuals.begin(), individuals.end(), list_bests.begin());

	   GetDistant distant(m_previous);
	   partial_sort(list_bests.begin(), list_bests.begin()+max, list_bests.end(), distant);
	   elem_dist = min_element(list_bests.begin(), list_bests.begin()+max, SortInd());
	   elem_id = (*elem_dist)->getId();
	}
    }

    return elem_id;
}

/**
void SelectXLS::getIndsToImprove(PopulationReal *pop, deque<tIndividualReal*> &subpop) {
    unsigned popsize = pop->size();
    tIndividualReal* ind;
    unsigned i, max;

    subpop.clear();

    // It select all individuals that have not been applied the LS and not improve
    for(i = 0, max=0; i < popsize; i++) {
	ind = pop->getInd(i);

	if (ind->getCount("non_improved")==0 && ind->isEval()) {
	    subpop.push_back(ind);
	    max++;
	}
    } 

}   

unsigned SelectXLS::selectIndToImprove(deque<tIndividualReal*> &individuals) {
    deque<tIndividualReal*>::iterator elem;
    deque<double> improvement;
    unsigned elem_id;

    assert(individuals.size()>0);
    
    if (m_previous.empty()) {
	elem = min_element(individuals.begin(), individuals.end(), isBetter);
	elem_id = (*elem)->getId();
    }
    else {
        elem = find_if(individuals.begin(), individuals.end(), isImproving);
	deque<tIndividualReal*>::iterator posi_max;

	if (elem != individuals.end()) {
	    elem_id = (*elem)->getId();
	}
	else {
	   vector<tIndividualReal*> list_bests(individuals.size());
	   vector<tIndividualReal*>::iterator elem_dist;
	   unsigned max = m_maxdists;

	   if (max > individuals.size()) {
	      max = individuals.size();
	   }


	   copy(individuals.begin(), individuals.end(), list_bests.begin());

	   GetDistant distant(m_previous);
	   partial_sort(list_bests.begin(), list_bests.begin()+max, list_bests.end(), distant);
	   elem_dist = min_element(list_bests.begin(), list_bests.begin()+max, SortInd());
	   elem_id = (*elem_dist)->getId();
	}
    }

    return elem_id;
}
*/
