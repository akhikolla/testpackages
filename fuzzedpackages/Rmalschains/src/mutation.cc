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

#include "mutation.h"
#include "problem.h"
#include "random.h"
#include <cassert>
#include <string>

using std::string;

using namespace realea;

Mutation::Mutation(IMutation *mut, double ratio) {
    assert(mut != NULL);
    m_mutation = mut;

    if (ratio < 0 || ratio > 1.0) {
	throw new string("Mutation::ratio is not a valide number");
    }

    m_ratio = ratio;
    m_random = NULL;
}

void Mutation::setRatio(double ratio) {
    if (ratio < 0 || ratio > 1.0) {
        throw new string("Mutation::ratio is not a valide number");
    }

    m_ratio = ratio;
}

bool Mutation::apply(tChromosomeReal &sol) {
    bool changed = false;
    unsigned size = sol.size();

    if (m_ratio > 0) {
    	if (!m_random) {
       	   throw ConfigException("Mutation::random");
    	}
 
	changed = (m_random->rand() <= m_ratio);
    }

    if (changed) {
        unsigned pos;

	do {
	    pos = m_random->randint(0, size-1);
	} while(!m_domain->canBeChanged(pos));

	sol[pos] = m_mutation->mutate(sol, pos);
    }

    return changed;
}

void Mutation::setRandom(Random *random) {
	m_random = random;
	m_mutation->setRandom(random);
}

void Mutation::setDomain(DomainRealPtr domain) {
	m_mutation->setDomain(domain);
	m_domain = domain;
}

tGen MutationBGA::mutate(tChromosomeReal &sol, unsigned pos) {
    const unsigned num=16;
    const double pai = 1.0/num;
    tReal rangi, min, max;
    unsigned i;
    double dif;
    double sum=0;

    if (!m_domain) {
	new ConfigException("MutationBGA::domain");
    }

    m_domain->getValues(pos, &min, &max);
    rangi = 0.1*(max-min);

    if (!m_random) {
	new ConfigException("MutationBGA::random");
    }

    for (i = 0, dif = 1; i < num; i++, dif/=2) {
	if (m_random->rand() < pai) {
	    sum += dif;
	}

    }

    tGen value = sol[pos];

    if (sum == 0) {
	return value;
    }

    // Obtain the sign
    if (m_random->rand() < 0.5) {
       value += rangi*sum;

       if (value > max && m_domain->isBound()) {
	  value = max;
       }
    }
    else {
	value -= rangi*sum;

	if (value < min && m_domain->isBound()) {
	    value = min;
	}
    
    }
 
    return(value);
}

Mutation::~Mutation(void) {
   if (m_mutation) {
      delete m_mutation;
   }
}
