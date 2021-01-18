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

#include "populationreal.h"
#include "random.h"
#include "initind.h"
#include "define.h"
#include "distance.h"
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <cassert>
#include <limits>

using std::sort;
using namespace realea;
using realea::internal::SimpleInitInd;
using realea::internal::UniformInitInd;

void PopulationReal::setObserver(PopulationObserver *observer) {
     m_observers.push_back(observer);
}

PopulationReal::PopulationReal(Random *random, unsigned int max, unsigned int pob): 
m_aditionalsize(max), m_size(pob), m_individuals(), 
m_knownbest(false), m_knownworst(false), m_observers(), m_random(random) {
    m_individuals.reserve(m_aditionalsize);
    m_initInd = new SimpleInitInd(m_random);
//    m_initInd = new UniformInitInd(m_random);
}

PopulationReal::~PopulationReal(void) {
    vector<tIndividualRealPtr>::iterator item;
    deque<PopulationObserver*>::iterator observer;


    for (item = m_individuals.begin(); item != m_individuals.end(); ++item) {
	delete (*item);
    }

    m_observers.clear();

    if (m_initInd) {
	delete m_initInd;
    }
}

template<class T>
void call_reset(T *elem) {
    elem->reset();
}

void PopulationReal::reset(DomainRealPtr domain, int posi) {
    unsigned size = m_individuals.size();
    tIndividualRealPtr old=NULL;

    m_knownbest = m_knownworst = false;

    if (posi >= 0) {
	old = m_individuals[posi];
	m_individuals[posi] = NULL;
    }

    for (unsigned i = 0; i < size; ++i) {
        if (m_individuals[i]) {
	   delete m_individuals[i];
	   m_individuals[i] = NULL;
	}
    }
    if (!m_individuals.empty()) {
       m_individuals.erase(m_individuals.begin(), m_individuals.end());
    }

    if (old != NULL) {
	m_initInd->resetBest(domain, old->sol(), m_size);
    }
    else {
	m_initInd->reset(domain, m_size);
    }

    for (unsigned i = 0; i < m_size; ++i) {
        tChromosomeReal crom;

	if (i != ( (unsigned) posi) ) {
	    m_initInd->createInd(domain, crom);
	    tIndividualReal *individual = getInstance(crom);
	    individual->setId(i);
	    m_individuals.push_back(individual);
	}
	else {
	    m_individuals.push_back(old);
	}
    }

    resetObservers();
}

tIndividualReal *PopulationReal::getInstance(tChromosomeReal &crom) {
    return new tIndividualReal(crom);
}

tIndividualReal *PopulationReal::getInstance(tChromosomeReal &crom, tFitness fitness) {
    return new tIndividualReal(crom, fitness);
}

bool isNull(tIndividualReal *ind) {
	return (ind==NULL);
}

bool PopulationReal::thereNull(void) {
    return false;
    vector<tIndividualReal*>::iterator posNull = find_if(m_individuals.begin(), m_individuals.end(), isNull);

    return (posNull != m_individuals.end());
}

void PopulationReal::updateObservers(void) {
     vector<tIndividualReal*>::iterator posind;
     deque<PopulationObserver*>::iterator observer;
     tIndividualRealPtr individual;
     unsigned newid, oldid;
     
     for (newid = 0, posind = m_individuals.begin(); posind != m_individuals.end(); ++posind) {
	for (observer = m_observers.begin(); observer != m_observers.end(); ++observer) {
	    individual = *posind; 
	    oldid = individual->getId();
	    newid++;
	    individual->setId(newid);
	    (*observer)->changeId(oldid, newid);
	}
     
     }
}

/**
 * Update the identifications of a specific individual
 *
 * @param id individual id to notify
 */
void PopulationReal::notifyObservers(unsigned id) {
     deque<PopulationObserver*>::iterator observer;
     
     for (observer = m_observers.begin(); observer != m_observers.end(); ++observer) {
	 (*observer)->notifyChange(id);
     }
}

/**
 * Update the identifications of individuals (it is called after removeWorses)
 * 
 * @see removeWorses sort
 */
void PopulationReal::resetObservers(void) {
     deque<PopulationObserver*>::iterator observer;
	for (observer = m_observers.begin(); observer != m_observers.end(); ++observer) {
	    (*observer)->reset();
     }
}



void PopulationReal::sort(void) {
    tIndividualReal::sort(m_individuals); 

    m_knownbest = m_knownworst = true;
    m_best = 0;
    m_worst = m_individuals.size()-1;

    if (thereNull()) {
	throw new string("PopulationReal::sort,there is null");
    }

}

bool isNotEval(tIndividualReal *ind) {
	return (!ind->isEval());
}

void PopulationReal::removePending(void) {
     if (m_individuals.empty())
	return;

    // All not evaluated are at the end
    vector<tIndividualReal*>::iterator beginNotEval = find_if(m_individuals.begin(), m_individuals.end(), isNotEval);

    // Check there is eval and not eval
    if ((beginNotEval != m_individuals.end())) {
	remove(beginNotEval-m_individuals.begin(), m_individuals.size());
     }
}

void PopulationReal::remove(unsigned begin, unsigned end) {
    assert(begin <= end);
    assert(end <= m_individuals.size());

    for (unsigned i = begin; i < end; i++) {
	delete m_individuals[i];
	m_individuals[i]=NULL;
    }

    vector<tIndividualReal*>::iterator base = m_individuals.begin();
    m_individuals.erase(base+begin, base+end);
}

void PopulationReal::removeWorses(void) {
    sort();
    removePending();
    unsigned size = m_individuals.size();

    if (!m_individuals.empty() && size > m_size) {
	remove(m_size, size);
    }

    updateObservers();
    m_worst = m_individuals.size()-1;
}

void PopulationReal::random(void) {
    unsigned int size=m_individuals.size();
    int ssize, pos;

    if (size == 0)
      return;

    m_knownbest = m_knownworst = false;

    ssize = size;
    int *sample = new int[size];
    initSample(sample, ssize);

    for (unsigned i = 0; i < size; ++i) {
	pos = m_random->getSample(sample, &ssize);
	swap(m_individuals[i], m_individuals[pos]);
    }

    delete[] sample;
}

tIndividualReal *PopulationReal::getInd(unsigned int pos) {
    assert(pos < m_size);
    return m_individuals[pos];
}

void PopulationReal::append(tChromosomeReal &sol, tFitness fitness) {
    tIndividualRealPtr ind = getInstance(sol, fitness);
    ind->setId(m_individuals.size());
    m_individuals.push_back(ind);
}

void PopulationReal::append(tIndividualReal *real) {
    if (m_individuals.size() == m_aditionalsize) {
	throw new runtime_error("maximum number of elems in population");
    }

    m_individuals.push_back(real);
    m_knownbest = m_knownworst = false;
}

void PopulationReal::change(unsigned int pos, const tChromosomeReal &sol, tFitness fitness) {
    assert(pos < m_size);
    m_individuals[pos]->change(sol, fitness);
    m_knownbest = m_knownworst = false;
}

unsigned PopulationReal::getBest(void) {
   assert(!m_individuals.empty());
   unsigned int ind;
   tIndividualReal *current, *best;
   int pos_best;

   if (m_knownbest) {
      return m_best;
   }

   pos_best = -1;

   best = NULL;
   unsigned size = m_individuals.size();
   
   for (ind = 0; ind < size; ++ind) {
	current = m_individuals[ind];

	if (!current->isEval()) {
	    continue;
	}

	if  (best == NULL) {
	    pos_best = ind;
	    best = current;
	}
	else if (current->isBetter(best)) {
	    best = current;
	    pos_best = ind;
	}
   }

   m_best = pos_best;
   m_knownbest = true;
   return pos_best;
}

class PopulationSort {
   public:    
    bool operator()(const unsigned &posa, const unsigned &posb) {
	if (!m_inds[posa]->isBetter(m_inds[posb])) 
	    return false;
	else 
	  return true;
    }
    PopulationSort(vector<tIndividualRealPtr> &v) : m_inds(v) {}
private:
    vector<tIndividualRealPtr> m_inds;
};

vector<unsigned> PopulationReal::getBests(unsigned num) {
    vector<unsigned> positions(m_size);
    vector<unsigned> bests(num);

    for (unsigned i = 0; i < m_individuals.size(); i++) {
	positions[i] = i;
    }
    
    partial_sort(positions.begin(), positions.begin()+num, positions.end(), 
		    PopulationSort(m_individuals));

    copy(positions.begin(), positions.begin()+num, bests.begin());
    return bests;
}

unsigned PopulationReal::getWorst(void) {
   assert(!m_individuals.empty());
   unsigned int ind;
   tIndividualReal *current,*worst;
   int pos_worst=-1;

    if (m_knownworst) {
	return m_worst;
    }

   worst = NULL;
   unsigned size = m_individuals.size();
   
   for (ind = 0; ind < size; ++ind) {
	current = m_individuals[ind];

	if (!current->isEval()) {
	    continue;
	}

	if (worst == NULL) {
	   pos_worst = ind;
	   worst = current;
	}
	else if (current->isWorse(worst)) {
	    worst = current;
	    pos_worst = ind;
	}
   }

   m_worst = pos_worst;
   m_knownworst = true;
   return pos_worst;
}


void PopulationReal::eval(IEvalInd *evalInd, unsigned neweval) {
    vector<tIndividualRealPtr>::iterator item;
    bool maxEval = false;

    for (item = m_individuals.begin(); item != m_individuals.end() && !maxEval; ++item) {
	if (!(*item)->isEval()) {
	   neweval -= 
	   evalInd->eval(*item);

	   if (neweval == 0) {
	      maxEval = true;
	   }
	}

    }
}

unsigned PopulationReal::size(void) {

    if (m_individuals.size() < m_size)
	return m_individuals.size();
    else 
        return m_size;
}

tFitness PopulationReal::getSecondBestFitness(void) {
    tFitness fitness, secondbestFitness=-1;
    unsigned best;
    int secondbest=-1;
    unsigned ind, popsize=size();
    
    best = getBest();

    for (ind = 0; ind < popsize; ++ind) {
       if (ind != best) {
	  fitness = m_individuals[ind]->perf();

	  if (secondbest < 0 || fitness < secondbestFitness) {
	     secondbest = ind;
	     secondbestFitness = fitness;
	  }
       }
    }

    assert(secondbest>=0);

    return secondbestFitness;
}

unsigned PopulationReal::ndim(void) {
    assert(m_individuals.size() > 0);
    return m_individuals[0]->sol().size();
}

void PopulationReal::replace(unsigned pos, tIndividualRealPtr newind) {
    tIndividualRealPtr old = m_individuals[pos];
    m_individuals[pos] = newind;
    m_individuals[pos]->setId(old->getId());
    delete old;

    notifyObservers(pos);
    
    // Update the best and worst positions
    if (m_knownbest) {
	if (m_best == pos) {
	    m_knownbest = false;
	}
	else if (newind->isBetter(m_individuals[m_best])) {
	    m_best = pos;
	}
    }

    if (m_knownworst) {
	if (m_worst == pos) {
	    m_knownworst = false;
	}
	else if (newind->isWorse(m_individuals[m_worst])) {
	   m_worst = pos;
	}
    }

}

void PopulationReal::replaceWithoutDeleting(unsigned pos, tIndividualRealPtr newind) {
    tIndividualRealPtr old = m_individuals[pos];

    m_individuals[pos] = newind;
    m_individuals[pos]->setId(old->getId());

    //delete old;

    notifyObservers(pos);

    // Update the best and worst positions
    if (m_knownbest) {
	if (m_best == pos) {
	    m_knownbest = false;
	}
	else if (newind->isBetter(m_individuals[m_best])) {
	    m_best = pos;
	}
    }

    if (m_knownworst) {
	if (m_worst == pos) {
	    m_knownworst = false;
	}
	else if (newind->isWorse(m_individuals[m_worst])) {
	   m_worst = pos;
	}
    }

}

tFitness PopulationReal::getMean(void) {
    tFitness sum=0;
    unsigned i;
    assert(!m_individuals.empty());

    for (i = 0; i < m_individuals.size() && m_individuals[i]->isEval(); ++i) {
	sum += m_individuals[i]->perf();
    }

    assert(i == m_individuals.size());
    return sum/m_individuals.size();
}

tFitness PopulationReal::getMedian(void) {
    vector<unsigned> positions(m_size);
    unsigned num, posi;

    for (unsigned i = 0; i < m_individuals.size(); i++) {
	positions[i] = i;
    }
 
    num  = m_individuals.size()/2;
    partial_sort(positions.begin(), positions.begin()+num, positions.end(), 
		    PopulationSort(m_individuals));

    posi = positions[num-1];
    return m_individuals[posi]->perf();
}

void PopulationReal::getPercentils(double *percen, unsigned num) {
    vector<unsigned> positions(m_size);
    unsigned i, posi;

    for (i = 0; i < m_individuals.size(); i++) {
	positions[i] = i;
    }
 
    partial_sort(positions.begin(), positions.end(), positions.end(), 
		    PopulationSort(m_individuals));

    percen[0] = m_individuals[getBest()]->perf();

    for (i = 1; i <= num; ++i) {
        posi = i*m_size/num;
	posi = positions[posi-1];
	percen[i] = m_individuals[posi]->perf();
    }

}

static PopulationObserver *g_observer=NULL;

static void assignd(vector<tIndividualReal*> &individuals, int dst, int orig) {
    if (individuals[dst] != NULL) {
    	delete individuals[dst];
    }

    individuals[dst] = individuals[orig];

    individuals[dst]->setId(dst);
    
    //TODO: Change when the LS Memory is not cleaned after reducing population
    //individuals[dst]->setId(dst);
    if (g_observer != NULL) {
        g_observer->changeId(dst, orig);
    }

    individuals[orig] = NULL;
}

void PopulationReal::reduceHalf(void) {
    int NP = m_individuals.size();
    int NewPopsize = (NP+1)/2;
    int i;

    if (!m_observers.empty()) {
	g_observer = m_observers.front();
    }

    for (i=0; i< NP/4; i++) {
        if (m_individuals[NP/4+i]->isBetter(m_individuals[i])) {
            assignd(m_individuals, i, NP/4+i);
        }
    }

    for (i=0; i< NP/4; i++) {
        if (m_individuals[3*NP/4+i]->isBetter(m_individuals[NP/2+i])) {
	    assignd(m_individuals, NP/2+i, 3*NP/4+i);
        }
    }

    for (i=0; i< NP/4; i++) {  // copy
 	assignd(m_individuals, NP/4+i, NP/2+i);
    }

    // Copy the last one if it was not divible by 2 (25 in the testing)
    if ((NP % 2)) {
	assignd(m_individuals, NewPopsize-1, NP-1);
    }

    for (i = NewPopsize; i < NP; i++) {
	if (m_individuals[i]!= NULL) {
	   delete m_individuals[i];
	   m_individuals[i] = NULL;
	}
    }
 
//    printf("Reduced: %d -> %d\n", NP, NewPopsize);

    m_individuals.erase(m_individuals.begin()+NewPopsize, m_individuals.end());

    m_knownbest = m_knownworst = false;
//    resetObservers();
}


double PopulationReal::getDiversity(void){



//    unsigned dim = m_individuals[0]->sol().size();
//
//    tChromosomeReal averagePoint(dim);
//
//    for(int i=0;i<dim;i++)
//    {
//        averagePoint[i]=0;
//        for(int j=0;j<m_size;j++)
//        {
//            averagePoint[i] += m_individuals[j]->sol()[i];
//        }
//        averagePoint[i] = averagePoint[i]/m_size;
//    }
//    double diversity = 0;
//    for(int j=0;j<m_size;j++)
//    {
//        diversity+=distreal(m_individuals[j]->sol(),averagePoint)*dim;
//    }
//
//    return diversity/m_size;


    double diversity=numeric_limits<double>::max();

    for(unsigned int i=0;i<m_size-1;i++)
    {
        for(unsigned int j=i+1;j<m_size;j++)
        {
            //printf("dist between %d and %d is %e\n", i,j,distBestSol(m_individuals[i],m_individuals[j]));
            diversity = min(distreal(m_individuals[i]->sol(),m_individuals[j]->sol()), diversity);
        }
    }
    return diversity*m_individuals[0]->sol().size();
}
