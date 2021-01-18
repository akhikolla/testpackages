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
#ifndef _PSO_H

#define _PSO_H 1

#include "iea.h"

namespace realea {

class ConfigPSO;

/**
 * @class Implementa el algoritmo PSO real
 */
class PSO : public ClassEAlgorithm {
public:
    PSO(Random *random);
    ~PSO(void);
    /**
     * Allow to set the Configuration for the individuals
     */
    void setConfigPSO(ConfigPSO *config);

    /**
     * By default the population size is 20 individuals
     */
    unsigned getDefaultPopsize(void);

    unsigned init(void);
    unsigned realApply(tChromosomeReal &sol, tFitness &fitness);
    void setPopsize(unsigned popsize);
private:
    ConfigPSO* m_config;
};

}
#endif
