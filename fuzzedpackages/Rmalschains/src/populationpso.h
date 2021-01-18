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
#ifndef _INDIVIDUALPSO_H

#define _INDIVIDUALPSO_H 1

#include "running.h"
#include "individual.h"
#include "domain.h"
#include "populationreal.h"
#include "random.h"

namespace realea {

class ConfigPSO {
   public:
       ConfigPSO(DomainRealPtr domain, double inergymin, double inergymax);
   public:
       double c1(void) { return m_c1; }
       double c2(void) { return m_c2; }
       double  x(void) { return m_x; }
       double w(double ratio);  
       double vmax(unsigned i) { return m_vmax[i]; }

   private: 
       double m_c1, m_c2, m_x;
       double m_wmax, m_wmin;
       tChromosomeReal m_vmax;
};

void setGlobalConfigPSO(ConfigPSO *config);
void delGlobalConfigPSO(void);

/**
 * This class store with the current best solution the current position and the velocity 
 */
class tIndividualPSO : public tIndividualReal {
public:
    /**
     * Create a new individual. 
     *
     * @param initialPos initial position of the particle
     *
     * @param random random generator (to create the velocity
     */
    tIndividualPSO(const tChromosomeReal &initialPos, Random *random);
    tIndividualPSO(const tChromosomeReal &initialPos, double fitness, Random *random); 
    /**
     * Set the fitness of the current best solution
     */
    void setCurrentFitness(double newfit);

    tChromosomeReal &current(void) {
	return m_current;
    }

    /**
     * Move the individual to next position
     *
     * @param best best individual to consider (attraction point)
     */
    void move(tChromosomeReal &best, double ratio);





private:
    /**
     * Init the PSO Parameters
     *
     * @param random random generator
     */
    void initPSOParams(Random *random); 

private:
    tChromosomeReal m_current;
    tChromosomeReal m_velocity;
};

/**
 * This class allow us to define a population of particles
 * 
 */
class PopulationPSO : public PopulationReal {
  public:
    PopulationPSO(Random *random,unsigned int max, unsigned int pob);
    void restart(DomainRealPtr domain);
    /**
     * Move all the individual and eval them using the evalInd object
     *
     * @param evalInd individual.
     * @param running stop criterion (to stop when it is achieved)
     */
    void move(IEvalInd *evalInd, Running *running);

   private:
    tIndividualReal* getInstance(tChromosomeReal &crom);
    tIndividualReal* getInstance(tChromosomeReal &crom, double fitness);
};

}

#endif
