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

#ifndef _IMA_H

#define _IMA_H 1

namespace realea {

/**
 * @class MemeticHart
 * @ingroup realea_common
 *
 * @brief Classical Memetic Algorithm (using the Hart's model)
 *
 * This class lets to apply the LS process to the current individual with a ratio 
 */
class MemeticHart : public ClassEAlgorithm {
public:
     
     MemeticHart(ClassEAlgorithm *alg, ILocalSearch *ls, unsigned maxitera, double ratio=0.0625) : m_alg(alg), m_ls(ls) {
     }

     /**
      * @return the default popsize
      */
    virtual unsigned getdefaultpopsize(void) {
	
    }

	
	
private:
	ClassEAlgorithm *m_alg;
	LocalSearch *m_ls;
	unsigned m_intensity;
	double m_ratio;
};

}

