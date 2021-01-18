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

#ifndef _INDIVIDUAL_H

#define _INDIVIDUAL_H 1

#include "real.h"
#include "define.h"
#include <map>
#include <deque>
#include <stdexcept>

using namespace std;

namespace realea {

class IndException: public runtime_error {
private:
	string m_name;

public:
	IndException(string name) : runtime_error(name) {
		m_name = name;
	}

public:
	virtual ~IndException(void) throw () { }
	virtual const char *what(void) {
		return m_name.c_str();
	}
};


typedef pair<string, unsigned> tCounter;
typedef deque< pair<string,unsigned> > tCounterList;



class tIndividualReal {

   public:
        tIndividualReal (const tChromosomeReal &com);
        tIndividualReal (const tChromosomeReal &com, tFitness fitness);


/**
 * Assign a change into the chromosome, indicating the fitness
 *
 * @param sol new solution
 * @param fitness new fitness
 * @param maintain true if maintain the identity (counters)
 */
	void change(const tChromosomeReal &sol, tFitness fitness); 

/**
 * Return the sol (to the crossover or the LS)
 * @return the chromosome
 */
	const tChromosomeReal &sol(void) {
	    return m_sol;
	}

/**
 * Return the performance
 * @return perf its performance measure (fitness)
 */
	tFitness perf(void);
	/**
	 * Increm a counter for that individual.
	 *
	 * This allow to count the number of events to the same individuals (selected,crossed,improvement by LS, ...)
	 *
	 * @param ident identificator of the counter
	 * @see getCount
	 */
	void incremCount(string ident);
	/**
	 * Return a counter value
	 * @param ident counter's identificator
	 * @return the count 
	 *
	 * @see incremCount
	 */
	unsigned int getCount (string ident);
	/**
	 * @return true is the individual has been evaluated
	 */
	bool isEval(void);

	/**
	 * Eval the individual
	 *
	 * @param IEval evaluation function
	 * @return fitness obtained
	 */
	tFitness eval(IEval *funeval);


	/**
	 * get the individual value
	 *
	 * @param pos (gen position, 0.. ndim-1)
	 * @return the current gen value
	 */
	tGen gen(unsigned pos);

	/**
	 * @return true if the current individual is better than the other
	 */
	bool isBetter(tIndividualReal *theother);

	/**
	 * @return true if the current individual is worse than the other
	 */
	bool isWorse(tIndividualReal *theother);

	/**
	 * Set the minimize criterion
	 */
	static void setMinimize(void);
	/**
	 * Set the minimize criterion
	 */
	static void setMaximize(void);

	/**
	 * Get the criterion
	 *
	 * @return true is setMinimize() has been set
	 */
	static bool isMinimize(void);

	/**
	 * Set the id
	 */
	void setId(unsigned id);

	unsigned getId(void);

	/**
	 * Destructor
	 */
	virtual ~tIndividualReal(void);

	/**
	 * Sort the vector of individuals
	 */
	static void sort(vector<tIndividualReal*> &individuals); 

protected:
	tChromosomeReal m_sol; /**< Solution represented by the individual */


   private:
	/**
	 * Set the additional value
	 *
	 * @param perf fitness value
	 * @param id identification
	 */
	void setPerf(tFitness perf);

	double m_id; /**< Identity id */

	tFitness m_perf; /**< Fitness */
	bool m_evaluated; /**< Store if the individual has been evaluated */
	deque< pair<string,unsigned> > pcounters; /**< Statistical counters */
	bool m_notid; /** Check the identity */
};

typedef tIndividualReal *tIndividualRealPtr;

/**
 * This class allows to evaluate an individual
 */
class IEvalInd : public IEval {
public:
    virtual unsigned eval(tIndividualReal *newind)=0;
    virtual tFitness eval(const tChromosomeReal &sol)=0;
};


}
#endif
