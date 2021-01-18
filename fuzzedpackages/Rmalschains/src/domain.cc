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

#include "domain.h"
#include <string>
#include <cassert>
#include <cstdio>

using namespace std;
using namespace realea;

DomainReal::DomainReal(unsigned int dim) : m_mins(dim), m_maxs(dim), m_dim(dim), m_isbound(true), m_search_ini(0), m_search_fin(dim-1) {
   m_check_dim = new bool[dim];

   fill(m_check_dim, m_check_dim+dim, true);
}

void DomainReal::checkGen(unsigned int gen) {
    if (gen >= m_dim) {
        char msg[100];
	sprintf(msg, "position %d is not valide for a gen", gen);
	throw new string(msg);
    }
 
}

void DomainReal::getValues(unsigned int gen, tReal *pmin, tReal *pmax, bool check) {
   if (check)
     checkGen(gen);

   *pmin = m_mins[gen];
   *pmax = m_maxs[gen];
}

void DomainReal::setValues(unsigned int gen, tReal min, tReal max,bool check) {
   if (check)
       checkGen(gen);

   assert(min <= max);

    m_mins[gen] = min;
    m_maxs[gen] = max;
}

tReal DomainReal::clip(unsigned int gen, tReal value, bool check) {
    tReal min, max;

    if (check)
       checkGen(gen);   

    if (!m_isbound) {
	return value;
    }

    min = m_mins[gen];
    max = m_maxs[gen];

    if (value < min)
    {
      return min;
    }
    else if (value > max)
    {
      return max;
    }
    else 
      return value;
}

bool DomainReal::canBeChanged(unsigned dim) {
    assert(dim < m_dim);

    return m_check_dim[dim];

    if (dim >= m_search_ini && dim <= m_search_fin) 
       return true;
    else
	return false;
}

void DomainReal::clip(tChromosomeReal &crom) {
    assert(crom.size() == m_dim);

    

    if (!m_isbound) {
	return;
    }

    for (unsigned i = 0; i < m_dim; ++i) {
	crom[i] = clip(i, crom[i], false);
    }
}

void DomainReal::clip(double *crom) {

    if (!m_isbound) {
	return;
    }

    for (unsigned i = 0; i < m_dim; ++i) {
	crom[i] = clip(i, crom[i], false);
    }
}



bool DomainReal::check(const tChromosomeReal &crom) {
    bool follow = true;
    assert(crom.size() == m_dim);

    for (unsigned i = 0; i < m_dim && follow; ++i) {
        if (crom[i] < m_mins[i])
	   follow = false;
	else if (crom[i] > m_maxs[i])
	   follow = false;
    }    
    return follow;
}

void DomainReal::setDomainCenter(tChromosomeReal center, double scale) {
    unsigned i;
    tReal min, max, range, value;

    for (i = 0; i < m_dim; i++) {
       getValues(i, &min, &max);
       value = center[i];
       range = (max-min)*scale/2;

       if (value - range > min) {
	  min = value - range;
       }
       if (value + range < max) {
	  max = value + range;
       }

       setValues(i, min, max);
    }

}


#include "random.h"

void DomainReal::getInitRandom(Random *random, tChromosomeReal &crom) {
    for (unsigned i = 0; i < m_dim; ++i) {
	crom.push_back(random->randreal(m_mins[i], m_maxs[i]));
    }
}

void DomainReal::getInit(Random *random, tChromosomeReal &crom) {
    getInitRandom(random, crom);
}

DomainReal::~DomainReal(void) {
    if (m_check_dim != NULL) {
	delete []m_check_dim;
    }
}

void DomainReal::setSearchDomain(bool *searchDim, unsigned dim) {
    assert(dim == m_dim);
    copy(searchDim, searchDim+dim, m_check_dim);
}

void DomainReal::getSearchDomain(bool *searchDim, unsigned dim) {
    assert(dim == m_dim);
    copy(m_check_dim, m_check_dim+dim, searchDim);
}
