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

#ifndef _CROSS_H

#define _CROSS_H 1

#include "domain.h"
#include "individual.h"
#include "random.h"
#include "signal.h"

#include "icross.h"

namespace realea {
/**
 * @class CrossBLX
 *
 * This class implements the Crossover operator BLX-\alpha
 *
 */
class CrossBLX : public ICrossBinary {
   public:
    /**
     * Constructor.
     *
     * @param alpha alpha value, its value sets the diversity. 0 => No diversity, 0.3 => Maintain diversity, 0.5-1 => Grow the diversity
     */
    CrossBLX(double alpha);

    /**
     * Destructor
     */
    virtual ~CrossBLX(void);
    
    void reset(void) {}
    /**
     * Create the operator
     */
    virtual void operator()(const tChromosomeReal &mom, tFitness fit_mom, const tChromosomeReal &dad, tFitness fit_dad, tChromosomeReal &child);


   private:
    double m_alpha;
};

/**
 * @class CrossPLX
 *
 * This class implements the Crossover operator PBLX-\alpha. It is like 
 * BLX but using the first parent (mother) as the center 
 * of the search. 
 * Also, it can be indicated the range of variable to be considered
 */
class CrossPBLX : public ICrossBinary {
   public:
    /**
     * Constructor.
     *
     * @param alpha alpha value, its value sets the diversity. 0 => No diversity, 0.3 => Maintain diversity, 0.5-1 => Grow the diversity
     *
     */
    CrossPBLX(double alpha);

    /**
     * Destructor
     */
    virtual ~CrossPBLX(void);
    
    void reset(void) {}
    /**
     * Create the operator
     */
    virtual void operator()(const tChromosomeReal &mom, tFitness fit_mom, const tChromosomeReal &dad, tFitness fit_dad, tChromosomeReal &child);

   private:
    double m_alpha;
    double m_dim_ini;
    double m_dim_fin;
};

/**
 * @class CrossDim
 *
 * This class implements the Crossover operator BLX-\alpha, only with a subsets of consecutive gens
 * Also, it can be indicated the range of variable to be considered
 */
class CrossDim : public ICrossBinary {
   public:
    /**
     * Constructor.
     *
     * @param alpha alpha value, its value sets the diversity. 0 => No diversity, 0.3 => Maintain diversity, 0.5-1 => Grow the diversity
     * @param p_r 
     *
     */
    CrossDim(double alpha, double p_r);

    /**
     * Destructor
     */
    virtual ~CrossDim(void);
    
    void reset(void) {}
    /**
     * Create the operator
     */
    virtual void operator()(const tChromosomeReal &mom, tFitness fit_mom, const tChromosomeReal &dad, tFitness fit_dad, tChromosomeReal &child);

   private:
    double m_alpha;
    double m_pr;
    unsigned ndim;
    double m_dim_ini;
    double m_dim_fin;
};



namespace internal {
/**
 * @class CrossBinary
 *
 * This class contain the real class that make the crossover, using a ICrossBinary that defines the operation
 */
class CrossBinary : public IReset {
public:
    CrossBinary(ICrossBinaryPtr cross);
    void operator()(tIndividualReal *mom, tIndividualReal *dad, tChromosomeReal &child);
    virtual void operator()(const tChromosomeReal &mom, tFitness fit_mom, const tChromosomeReal &dad, tFitness fit_dad, tChromosomeReal &child);
    virtual ~CrossBinary(void);
private:
    ICrossBinaryPtr m_cross;
};

typedef CrossBinary* CrossBinaryPtr;
}

}

#endif
