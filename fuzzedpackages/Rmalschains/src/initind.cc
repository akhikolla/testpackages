#include "initind.h"
#include <cmath>
#include <cassert>

#define EPSILON 1E-18

using namespace realea;
using namespace realea::internal;

void SimpleInitInd::createInd(DomainRealPtr domain, tChromosomeReal &crom) {
    domain->getInit(m_random, crom);
}

void SimpleInitInd::reset(DomainRealPtr domain, unsigned count) {
}

void SimpleInitInd::resetBest(DomainRealPtr domain,const tChromosomeReal &sol, unsigned count) {
}

ElemDimInit::ElemDimInit(tGen min, tGen max, unsigned count, unsigned intervals) : m_min(min), m_size(intervals), m_interval(intervals) {
    m_range = (max-min)/m_size;

    for (unsigned i = 0; i < m_size; i++) {
	ElemRangeInit elem;
	elem.interval = i;
	elem.count = count;
	m_interval[i] = elem;
    }

}

void ElemDimInit::reduce(tGen value) {
    unsigned pos_interval;

    pos_interval = (unsigned) floor((value - m_min)/m_range+EPSILON);

    if (pos_interval >= m_interval.size() ) {
       pos_interval = m_interval.size()-1;
    }
    m_interval[pos_interval].count--;
    assert(m_interval[pos_interval].count > 0);
}

tGen ElemDimInit::random(Random *random) {
    unsigned pos_interval = random->randint(0, m_size-1);
    unsigned interval = m_interval[pos_interval].interval;
    assert(m_interval[pos_interval].count > 0);

    tGen value = m_min + interval*m_range+random->randreal(0, m_range);


    m_interval[pos_interval].count--;
    
    if (m_interval[pos_interval].count == 0) {
	m_interval[pos_interval] = m_interval[m_size-1];
	m_size -= 1;
    }

    return value;
}


void UniformInitInd::reset(DomainRealPtr domain, unsigned count) {
    unsigned ndim = domain->getDimension(); 
    tGen min, max;
 
    m_interval_dim.clear();
    unsigned intervals = count/10;

    for (unsigned i = 0; i < ndim; ++i) {
	domain->getValues(i, &min, &max);
        ElemDimInit elem(min, max, 10, intervals);
	m_interval_dim.push_back(elem);
    }

}

void UniformInitInd::resetBest(DomainRealPtr domain, const tChromosomeReal &best, unsigned count) {
    deque<ElemDimInit>::iterator elemDim;
    unsigned i;
    reset(domain, count);

    for (i = 0, elemDim = m_interval_dim.begin(); elemDim != m_interval_dim.end(); elemDim++, i++) {
	elemDim->reduce(best[i]);
    }
}

void UniformInitInd::createInd(DomainRealPtr domain, tChromosomeReal &crom) {
    unsigned i;   
    deque<ElemDimInit>::iterator elemDim;

    for (i = 0, elemDim = m_interval_dim.begin(); elemDim != m_interval_dim.end(); elemDim++, i++) {
	crom.push_back(elemDim->random(m_random));
    }
}

UniformInitInd::~UniformInitInd(void) {
}
