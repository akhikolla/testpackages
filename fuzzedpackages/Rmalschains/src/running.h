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

#ifndef _RUNNABLE_H
#define _RUNNABLE_H 1

#include "signal.h"
#include "define.h"
#include <ctime>
#include <stdexcept>

using namespace std;

/**
 * @class OptimeCriterion
 * @ingroup realea_common
 *
 * @brief Let specify when a class is optimum.
 *
 * This class' target is easier to specify when the optimum is achieved
 * 
 */

#define MINIMIZE minimize
#define MAXIMIZE maximize

namespace realea {


struct OptimeCriterion { 
public:

   /**
    * Constructor.
    *
    * @param optime optimum value
    * @param dif threshold value
    */
   OptimeCriterion(double optime, double dif);
   /**
    * Set the threshold required to identify two solution are equivalent
    * @param dif diference value
    */
   void setThreshold(double dif);

   /**
    * @return the threshold
    */
   inline double getThreshold(void);

    /**
     * Check if the fitness is too close to optimum
     *
     * @param fitness fitness of current solution
     * @return true if abs(fitness - fitness_optimum) < threshold
     */
    bool isOptime(double fitness); 

	/**
	 * Set the minimize criterion
	 */
    void setMinimize(void) {
	m_minimize = true;
    }
	/**
	 * Set the minimize criterion
	 */
    void setMaximize(void) {
        m_minimize = false;
    }

    bool minimize(void) {
	return m_minimize;
    }

    bool maximize(void) {
	return !m_minimize;
    }


   /**
    * return if fitness value1 is better than value2
    *
    * @param value1
    * @param value2
    * @return true if value is better than value2
    */
   bool isBetter(double value1, double value2) {
       if (m_minimize) {
          return (value1 < value2);
       }
       else {
          return (value1 > value2);
	}

}

   double getOptime(void) {
	return m_optime;
   }

private:
   double m_optime;
   double m_threshold; /**< Threshold value */
   double m_minimize;
};

class RunningException : public runtime_error {
public:
    RunningException(string msg) : runtime_error(msg), m_msg(msg) {
    }

    const char *what(void) {
	return m_msg.c_str();
    }

    ~RunningException(void) throw () {}
private:
	string m_msg;
};

/**
 * @class Running
 *
 * @brief Allow to control the stopping criterion
 *
 * It detect if the algorithm must stop, or because the maximum evaluation number is achieved 
 * or because the difference between the current solution and optimum is < threshold
 */
class Running : public IReset, public IFinish {
public:
    /**
     * Increment the evaluation number 
     */
    void increm(void); 

 
    void notify(double fitness);

    /**
     * Reset the counting (for restart the experimentation)
     */
    void reset(void); 

    /**
     * Constructor
     * @param isOptime optimum criterion 
     */
    Running(OptimeCriterion *isOptime);
   /**
     * Set the maximum evaluation  
     *
     * @param maxEval new maximum evaluation number
     */
    void setMaxEval(unsigned int maxEval); 

    virtual ~Running(void);

    /**
     * Return a new subrunning (to let algorithm that uses anothers)
     *
     * @param submaxeval maxeval of the new running returning
     *
     * @return a new running with the maxeval indicated
     */
    Running*getSubRunning(unsigned submaxeval);

    /**
     * @return the maximum evaluation number
     */
    unsigned int maxEval(void); 
    /**
     * @return true if the algorithm must stop (following the optimum criterion of the maximum evaluations
     * number)
     */
    bool isFinish(void);

    /**
     * Set the ratio of maximum evaluation done
     *
     * @return a real value between 0 and 1, the ratio of number evaluation currency done 
     */
    double ratio();

    /**
     * Check if the fitness is too close to optimum
     *
     * @param fitness fitness of current solution
     * @return true if abs(fitness - fitness_optimum) < threshold
     */
    bool isOptime(double fitness); 

   /**
    * Set the threshold required to identify two solution are equivalent
    * @param dif diference value
    */
   void setThreshold(double dif);

   void setMaxTime(unsigned seconds);

    /**
     * @return the current evaluation number
     */
    unsigned int numEval(void);

    /**
     * @return the threshold used
     */
    double getThreshold(void); 

    /**
     * @return the maximum evaluation
     */
    unsigned pending(void) {
	return (m_maxeval-m_neval);
    }

    /**
     * It is notify when a new chromosome is evaluated
     *
     * @param fit new fitness obtained.
     */
    void notifyEval(double fit);
    /**
     * @return true if the optime has been achieved
     */
    bool hasFoundOptime(void) {
	return m_optimized;
    }

    bool isBetter(double fit1, double fit2) {
	return m_checkOptime->isBetter(fit1,fit2);
    }

private:
    unsigned int m_neval;  /**< Current evaluation number */

    unsigned int m_maxeval;  /**< Maximum evaluation number */

    int m_maxmsecs; /**< Maximum time number */
    clock_t m_timeInit; 

    OptimeCriterion *m_checkOptime;    /**< Optimum criterion applied */

    bool    m_optimized;  /**< Stores if the optimum has been achieved */

    Running *m_parent; /**< Stores a new Running to assign */
    list<Running*> m_children;

    double m_best;  /** @var Best fitness obtained */
};

typedef Running* RunningPtr;
}

#endif
