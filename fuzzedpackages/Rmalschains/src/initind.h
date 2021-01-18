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
#ifndef _INIT_IND_H 

#define _INIT_IND_H 1

#include "real.h"
#include "domain.h"
#include <deque>

using std::deque;

namespace realea {
/**
 * @class 
 * Interfaz that allow to create new individuals
 */
class InitIndividual {
public:
    /**
     * Reset the class (it is called after each EA run
     */
    virtual void reset(DomainRealPtr domain, unsigned count)=0;
    /**
     * Reset the class (it is called when the EA is reset maintaining the best individual)
     * @param best the best individual
     */
    virtual void resetBest(DomainRealPtr domain,const tChromosomeReal &best, unsigned count)=0;
    /**
     * Create a new individual
     * @param domain Domain of the problem.
     * @param crom new chromosome
     */ 
    virtual void createInd(DomainRealPtr domain, tChromosomeReal &crom)=0;
    virtual ~InitIndividual(void) {}
};

namespace internal {

class SimpleInitInd : public InitIndividual {
public:
    SimpleInitInd(Random *random) {
	m_random = random;
    }
    void reset(DomainRealPtr domain, unsigned count);
    void resetBest(DomainRealPtr domain,const tChromosomeReal &best, unsigned count);
    void createInd(DomainRealPtr domain, tChromosomeReal &crom); 
    ~SimpleInitInd(void) {}
private:
    Random *m_random;
};

struct ElemRangeInit {
    unsigned interval;
    unsigned count;
};

/**
 * @class
 *
 * This class allow to generate the random value for one dimension assuring that for there are distributed
 * equally for each interval
 *
 */
class ElemDimInit {
private:
    tGen m_min;
    unsigned m_size;
    vector<ElemRangeInit> m_interval;
    double m_range;
public:
    /**
     * Constructor
     * @param min minimum value
     * @param max maximum value
     * @param count maximum number of values for each interval
     * @param intervals interval number
     */
    ElemDimInit(tGen min, tGen max, unsigned count, unsigned intervals); 
    /**
     * Reduce the count with the existing value
     */
    void reduce(tGen value);
    /**
     * Generate a new random value considering the different intervals
     * @parma random random generator
     */
    tGen random(Random *random);
};

class UniformInitInd : public InitIndividual {
public:
    UniformInitInd(Random *random) : m_random(random), m_interval_dim() {
	m_random = random;
    }

    void reset(DomainRealPtr domain, unsigned count);
    void resetBest(DomainRealPtr domain,const tChromosomeReal &best, unsigned count);
    void createInd(DomainRealPtr domain, tChromosomeReal &crom); 
    ~UniformInitInd(void);
private:
    Random *m_random;
    deque<ElemDimInit> m_interval_dim;
};

}

}

#endif
