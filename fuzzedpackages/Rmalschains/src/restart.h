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
#ifndef _RESTART_H 

#define _RESTART_H 1

#include "real.h"
#include "define.h"
#include "problem.h"
#include "populationreal.h"

namespace realea {

/**
 * @class RestartStrategy Set the restart strategy over a population
 * 
 * It 
 */
class RestartStrategy {
public:
    virtual void apply(PopulationReal *pop_alg, Problem *problem, IEvalInd *eval)=0;
    virtual ~RestartStrategy() {}
};

/**
 * @class RestartBest. Restart all the population individuals except the best one
 */
class RestartBest: public RestartStrategy {
   public:
    RestartBest(void) {}
    void apply(PopulationReal *pop_alg, Problem *problem, IEvalInd *eval);
};

/**
 * @class RestartReduce. Restart all the population individuals except the best one, and set the
 * domain space
 */
class RestartReduce: public RestartBest {
   public:
    RestartReduce(double scale);
    void apply(PopulationReal *pop_alg, Problem *problem, IEvalInd *eval);
   private:
    double m_scale;
};


}
#endif
