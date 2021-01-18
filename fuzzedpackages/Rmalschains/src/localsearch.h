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

#ifndef _LOCALSEARCH_H

#define _LOCALSEARCH_H 1

#include "ilocalsearch.h"

/**
 * @ingroup realea::common::internal
 */
namespace realea {

namespace internal {

/**
 * @class NewIndividualLocalSearchManager
 *
 * It is the responsable of applied the Local Search when it is created a new solution
 *
 * This is responsable of calling to LS method when it is needed.
 */
class NewIndividualLocalSearchManager: public IReset {
   public:
    /**
     * Constructor
     */
     NewIndividualLocalSearchManager(ILocalSearch *ls);

     /**
      * Restart the local search process
      */
     void reset(void);

     /**
      * This method is called by each new individual, and applies the LS if it is needed.
      *
      * This method is called by each new individual generated. If must decide if apply the
      * LS or not, in the last case it called the method m_localsearch->apply().
      *
      *
      * @param sol new Chromosome to improve, it could ber modified
      * @param pfitness fitness of the solution, it will be updated if sol is changed
      *
      * @param options options to LS method (by default it is NULL).
      *
      * @return bool true if the current sol has been changed
      *
      * @see ILocalSearch::apply
      */
     virtual bool applyNewSol(tChromosomeReal &sol, tFitness *pfitness, ILSParameters *params=NULL) = 0; 

     /**
      * Destructor. Only removes the LocalSearch object.
      */
     virtual ~NewIndividualLocalSearchManager(void);
   protected:
	ILocalSearch *m_ls;
};

/**
 * @class LocalSearchClasicalManager
 *
 * Apply the LS method only to a ratio of new individuals with a const intensity
 */
class RatioLocalSearchManager : public NewIndividualLocalSearchManager {
   public:
      /**
       * Constructor
       *
       * @param ls local search method
       * @param intensity intensity of LS
       * @param ratio ratio of new individuals to be improved
       */
      RatioLocalSearchManager(ILocalSearch *ls, unsigned intensity, double porcen=0.0625);

     /**
      * @return the current intensity
      * @see setIntensity
      */
     void setRandom(Random *random);



      /**
       * Applied the LS only to ratio of new individuals
       *
       * @param sol new Chromosome to improve, it could ber modified
       * @param pfitness fitness of the solution, it will be updated if sol is changed
       *
       * @param options options to LS method (by default it is NULL).
       *
       * @return bool true if the current sol has been changed
       *
       * @see NewIndividualLocalSearchManager::applyNewSol
       */
      bool applyNewSol(tChromosomeReal &sol, tFitness *pfitness, ILSParameters *params=NULL);


   private:

     /**
      * @return the current intensity
      * @see setIntensity
      */
     unsigned getIntensity(void);


	Random *m_random; /**< Random generator */
        double m_ratio; /**< Ratio of new individuals to be improved */

 	unsigned m_intensity; /** Intensity of each LS */
	bool m_configured; 
};

}}

#endif
