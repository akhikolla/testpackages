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

#include "individual.h"
#include <algorithm>
#include "debug.h"
#include <cassert>
#include <cstdio>

using namespace realea;

static bool m_criterion = false;
static bool m_minimize;


/**
 * Sort the individuals (not evaluated at the end)
 */
class SortIndMin {
    public:
	bool operator()(tIndividualReal *a, tIndividualReal *b) {
	    if (a->isEval() && b->isEval())
		return (a->perf() < b->perf());
	    else if (a->isEval()) 
		return true;
	    else
		return false;
	}
};

class SortIndMax {
    public:
	bool operator()(tIndividualReal *a, tIndividualReal *b) {
	    if (a->isEval() && b->isEval())
		return (a->perf() > b->perf());
	    else if (a->isEval())
		return true;
	    else 
		return false;
	}
};



/**
 * Constructor de un individuo
 */
tIndividualReal::tIndividualReal(const tChromosomeReal &com) : m_sol(com) , m_evaluated(false), pcounters(),m_notid(true) {
   
}

tIndividualReal::tIndividualReal(const tChromosomeReal &com, tFitness fitness) : m_sol(com) , m_evaluated(true), pcounters(),m_notid(true) {
    m_perf  = fitness;
}


tIndividualReal::~tIndividualReal(void) {
    pcounters.clear();
}

void tIndividualReal::setMinimize(void) {
	m_minimize = true;
	m_criterion = true;
}


void tIndividualReal::setMaximize(void) {
	m_minimize = false;
	m_criterion = true;
}

bool tIndividualReal::isMinimize(void) {
	return m_minimize;
}


void tIndividualReal::change(const tChromosomeReal &sol, tFitness fitness) {
    m_evaluated = true;
    m_sol = sol;
    m_perf = fitness;
}


tFitness tIndividualReal::perf(void) {
    if (!m_evaluated) 
      throw new IndException("Performance measure has not been obtained");

    return m_perf;
}

struct FindKey : public unary_function<tCounter,bool> {
   string m_key;

    bool operator()(const tCounter &value) {
	return (value.first == m_key);
    }
    
};

void tIndividualReal::incremCount(string ident) {
    tCounterList::iterator posi;
    FindKey find_key;
    find_key.m_key = ident;

    posi = find_if(pcounters.begin(), pcounters.end(), find_key);

    if (posi==pcounters.end()) {
       tCounter count;
       count.first = ident;
       count.second = 1;

//       printf("individuo %u debe incrementar %s\n", getId(), ident.c_str());
       pcounters.push_back(count);
    }
    else {
	posi->second ++;
    }

}


unsigned int tIndividualReal::getCount(string ident) {
    tCounterList::iterator posi;
    FindKey find_key;
    find_key.m_key = ident;

    posi = find_if(pcounters.begin(), pcounters.end(), find_key);

   if (posi == pcounters.end())
      return 0;
   else 
      return posi->second;
}

void tIndividualReal::setPerf(tFitness perf) {
    if (m_evaluated) {
	throw new string("individual has been evaluated previously");
    }
    m_perf = perf;
    m_evaluated = true;
}

bool tIndividualReal::isEval(void) {
    return m_evaluated;
}

bool tIndividualReal::isBetter(tIndividualReal *theother) {
    if (!m_criterion) 
      throw new IndException("Criterion (Maximize/Maximize) has not been set");

    if (m_minimize) {
	return (this->perf() < theother->perf());
    }	
    else 
	return (this->perf() > theother->perf());
}

bool tIndividualReal::isWorse(tIndividualReal *theother) {
    if (!m_criterion) 
      throw new IndException("Criterion (Maximize/Maximize) has not been set");

    if (m_minimize) {
	return (this->perf() > theother->perf());
    }	
    else 
	return (this->perf() < theother->perf());

}

/**
 * Set the minimize criterion
 */
void setMinimize(void);
/**
 * Set the minimize criterion
 */
	void setMaximize(void);


tGen tIndividualReal::gen(unsigned n) {
    if (!(n < m_sol.size())) {
	print_error("Size: %u\tn: %u\n", m_sol.size(), n);
    }

    assert(n < m_sol.size());
    return m_sol[n];
}

tFitness tIndividualReal::eval(IEval *funeval) {
    if (m_evaluated) {
	return m_perf;
    }

    tFitness fitness = funeval->eval(this->sol());
    setPerf(fitness);
    return fitness;
}

void tIndividualReal::setId(unsigned id) {
    m_id = id;
    m_notid = false;
}

unsigned tIndividualReal::getId(void) {
   if (m_notid) {
      throw string("IndividualReal::Id");
   }
   return m_id;
}

void tIndividualReal::sort(vector<tIndividualReal*> &individuals) {
    if (tIndividualReal::isMinimize()) 
    	std::sort(individuals.begin(), individuals.end(), SortIndMin());
    else 
    	std::sort(individuals.begin(), individuals.end(), SortIndMax());
}
