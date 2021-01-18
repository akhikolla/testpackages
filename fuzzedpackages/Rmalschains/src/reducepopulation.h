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
#ifndef _REDUCEPOPULATION_H 

#define _REDUCEPOPULATION_H 1

#include "populationreal.h"
#include "running.h"


namespace realea {
/**
 * @class Implement the population reduction used in several DEs, proposed in jDEdynNP-FCR @see Jde
 *
 * It defines the required methods for population size, it follows the strategy pattern.
 */
class PopulationReductionStrategy {

public:
    /**
     * Constructor (disable by default)
     */
    PopulationReductionStrategy(void);

    void setNumReductions(int numReductions);

    /**
     * Set the number of reductions. 
     *
     * @param m_running Running criterion (to get the current evaluations number, and the ºmaximum FEs).
     */
    void config(Running *running);

    /*
     * reduce to half the popsize of pop if it is decided. 
     *
     * It reduces the population each maxFEs/(Reductions+1), to make at the end only numReductions.
     *
     * @param pop population to check and update.
     *
     * @return true if the population is changed.
     *
     */
    bool updatePopulationSize(PopulationReal *pop);

private: 
      Running *m_running;
      // Number of popsize reductions 
      int m_numReductions;
      int m_maxReductions;

      int m_countdown;

      /// obtained from m_running
      int m_initEvals; 
      int m_maxFES;
};

}

#endif
