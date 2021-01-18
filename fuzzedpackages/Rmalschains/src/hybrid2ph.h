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
#ifndef _MA2PH_H

#define _MA2PH_H 1

#include "hybrid.h"

namespace realea {

/**
 * @class MA2PH
 *
 * @brief hybrid algorithm using a two phase scheme (really it is not a MA)
 *
 * It is applied first the AE, and after the LS over the best one
 */
class Hybrid2Ph : public Hybrid {
public:
	/**
	 * Constructor
	 *
	 * @param alg EA to applied
	 * @param ls Local Search to apply
	 * @param intensity LS intensity (in number evaluations)
	 * @param ratio ratio of LS applications
	 */
	Hybrid2Ph(IEAlgorithm *alg, ILocalSearch *ls) : Hybrid(alg,ls) {
		m_effort = -1;
	}

	void setRunning(Running *running); 

	/**
	 * @param ratio. Set the global ratio invested into the Local Search
	 *
	 * @param ratio global ls/total ratio
	 * (change the intensity, because only one LS improvement is made)
	 */
	void setEffortRatio(double ratio);

	/**
	 * Set first the EA and after the LS (only one time, to best one)
	 */
	unsigned realApply(tChromosomeReal &sol, tFitness &fitness);

	void setMaxEval(unsigned int maxeval);

	unsigned init(void);

private:
	double m_effort; 
};

}

#endif
