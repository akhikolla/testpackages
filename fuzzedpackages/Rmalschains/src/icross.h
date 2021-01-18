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

#ifndef _ICROSS_H

#define _ICROSS_H 1

#include "domain.h"
#include "individual.h"
#include "random.h"
#include "running.h"
#include "signal.h"

namespace realea {
/**
 * @class ICross
 * @ingroup realea_common
 *
 * This class allow to create new individuals from two individuals
 *
 * To create new crossover operator this class must be inherited and defined 
 * the method operator.
 *
 */
class ICrossBinary : public IReset {
public:
    /**
     * Allow to select the domain related to the problem
     *
     * @param domain domain of the different values
     */
    void setDomain(DomainRealPtr domain) {
	m_domain = domain;
    }

    /**
     * Allow to set the random generator
     *
     * @param random 
     */
    void setRandom(RandomPtr random) {
	m_random = random;
    }

    /**
     * Allow to set the running generator
     *
     * @param running
     */
    void setRunning(RunningPtr running) {
	m_running = running;
    }

    /**
     * Generate the new chromosome child from two existing chromosomes.
     * (when this method is called, variables m_random and m_domain have been properly initialised).
     *
     * @param mom chromosome 'mother'
     * @param dad chromosome 'father'
     *
     * @param child the new chromosome, output
     *
     */
    virtual void operator()(const tChromosomeReal &mom, tFitness fit_mom, const tChromosomeReal &dad, tFitness fit_dad, tChromosomeReal &child)=0;


protected:
    DomainRealPtr m_domain; /**< Domain of variables */
    RandomPtr m_random;    /**< Random generator */
    RunningPtr m_running;  /**< Running evaluator, to allow it adapt their values */
};

typedef ICrossBinary* ICrossBinaryPtr;
}

#endif
