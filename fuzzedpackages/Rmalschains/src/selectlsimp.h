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
 *
 */

#ifndef _SELECT_LS_IMP_H

#define _SELECT_LS_IMP_H 1

#include "selectls.h"
#include <deque>

using std::deque;

namespace realea {

class SelectBestToImprove :  public SelectImprovementLS {
   public:
	/**
	 * Return the individual with the best fitness
	 *
	 * @param pop population to search
	 *
	 * @param subpop an deque of individuals with the elements to improve
	 */
	void getIndsToImprove(PopulationReal *pop, deque<tIndividualReal*> &subpop);

	/**
	 * Select all the individuals of the population that has not been selected or the last time it was improved
	 *
	 * @param individuals individuals to consider
	 *
	 * @return identification of the chosen individual 
	 *
	 */
	unsigned selectIndToImprove(deque<tIndividualReal*> &individuals); 
	~SelectBestToImprove(void) {}
};

class SelectWithDiversityToImprove:  public SelectImprovementLS {
   public:
       SelectWithDiversityToImprove(void) : m_previous() {}

	/**
	 * Select all the individuals of the population that has not been selected or the last time it was improved, 
	 * and store in the list the previously selected
	 * @param individuals individuals to consider
	 *
	 * @return identification of the chosen individual 
	 *  
	 */
	void getIndsToImprove(PopulationReal *pop, deque<tIndividualReal*> &subpop);

   protected:
	deque<tIndividualReal*> m_previous;
};


class SelectDiverseToImprove:  public SelectWithDiversityToImprove {
   public:
       SelectDiverseToImprove(void): SelectWithDiversityToImprove() {}

	/**
	 * If no individual has been improvement return the individual with the best fitness, 
	 * othercase it select the most distance to the previously selected
	 *
	 * @param pop population to search
	 *
	 * @param subpop an deque of individuals with the elements to improve
	 */

	unsigned selectIndToImprove(deque<tIndividualReal*> &individuals); 

};

class SelectDistantBestToImprove:  public SelectWithDiversityToImprove {
   public:
       SelectDistantBestToImprove(unsigned maxbest) : SelectWithDiversityToImprove(), m_maxbests(maxbest) {}
	/**
	 * Select all the individuals of the population that has not been selected or the last time it was improved, 
	 * and store in the list the previously selected

	 * @param individuals individuals to consider
	 *
	 * @return identification of the chosen individual 
	 *  
	 */
	unsigned selectIndToImprove(deque<tIndividualReal*> &individuals); 
   private:
	unsigned m_maxbests;
};

class SelectBestDistantToImprove:  public SelectWithDiversityToImprove {
   public:
       SelectBestDistantToImprove(unsigned maxdist) : SelectWithDiversityToImprove(), m_maxdists(maxdist) {}
	/**
	 * Return the individual with best fitness from the m_maxdists more distant solutions
	 *
	 * @param individuals individuals to consider
	 *
	 * @return identification of the chosen individual 
	 *  
	 */
	unsigned selectIndToImprove(deque<tIndividualReal*> &individuals); 
   private:
	unsigned m_maxdists;
};

/*
class SelectXLS:  public SelectImprovementLS {
   public:
       SelectXLS(void) : m_previous() {}

	void getIndsToImprove(PopulationReal *pop, deque<tIndividualReal*> &subpop);
	unsigned selectIndToImprove(deque<tIndividualReal*> &individuals); 

   protected:
	deque<tIndividualReal*> m_previous;
};
*/

}

#endif
