/*--------------------------------------------------------------------*/
/*     Copyright (C) 2013-2013  Serge Iovleff, Quentin Grimonprez

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : quentin.grimonprez@inria.fr
*/

/*
 * Project:  MPAGenomics::
 * created on: 4 nov. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file LogisticLasso.h
 *  @brief In this file .
 **/


#ifndef LOGISTICLASSO_H_
#define LOGISTICLASSO_H_


#include "PenalizedModels.h"
#include "LassoPenalty.h"
#include "LogisticLassoSolver.h"

namespace HD
{
  class LogisticLasso;

  template<>
  struct ModelTraits<LogisticLasso>
  {
    typedef LogisticLassoSolver Solver;
    typedef LassoPenalty Penalty;
    typedef LassoMultiplicator Multiplicator;
    typedef STK::CG<LassoMultiplicator,STK::VectorX, InitFunctor> CG;
  };

  /**
   * Class Lasso derived from @c PenalizedModels.
   * This class constructs a Lasso Model to be solved by an @c EM algorithm.
   */
  class LogisticLasso : public PenalizedModels<LogisticLasso>
  {
    public:
      typedef STK::CG<LassoMultiplicator,STK::VectorX, InitFunctor> CGSolver;
      /**
       * Constructor
       * @param p_x pointer to the data
       * @param p_y pointer to the response
       * @param lambda value of parameter associated to the l1 penalty
       * @param threshold threshold for setting coefficient to 0
       * @param epsCG epsilon for CG convergence
       */
      LogisticLasso( STK::ArrayXX const* p_x, STK::VectorX const* p_y
                   , STK::Real lambda, STK::Real threshold, STK::Real epsCG)
                   : PenalizedModels<LogisticLasso>(p_x, p_y)
      {
#ifdef HD_DEBUG
        std::cout << "Creating LogisticLasso.\nLambda =" << lambda << std::endl;
        std::cout << "threshold =" << threshold << std::endl;
        std::cout << "epsCG =" << epsCG << std::endl;
#endif
        // creation lasso penalty
        p_penalty_ = new LassoPenalty(lambda);
        //create solver for lasso and add the penalty
        p_solver_ = new LogisticLassoSolver(p_x_, p_y_, &beta_, threshold, epsCG, p_penalty_);
#ifdef HD_DEBUG
        std::cout << "LogisticLasso Initialized. Likelihood =" << p_solver_->computeLlc() << std::endl;
#endif
      }
      /** destructor*/
      inline virtual ~LogisticLasso() {}
      /** set the lasso regularization parameter
       *  @param lambda
       */
      inline void setLambda(STK::Real lambda) { p_penalty_->setLambda(lambda);}
      /** set the threshold for the shrinkage
       *  @param threshold
       */
      inline void setThreshold(STK::Real threshold) { p_solver_->setThreshold(threshold);}
      /** set the epsilon for the CG
       *  @param epsilon epsilon for the convergence of CG
       */
      inline void setCGEps( STK::Real eps) { p_solver_->setCGEps(eps);}
      /** initialization of the class with a new beta0 */
      void initializeBeta()
      {
        p_solver_->disruptsBeta();
        p_solver_->initializeSolver();
#ifdef HD_DEBUG
        std::cout << "LogisitcLasso::initializeBeta done. Likelihood =" << p_solver_->computeLlc() << std::endl;
#endif
      }
      /** initialize the containers of the solver*/
      void initializeModel()
      {
        //set the parameter of the solver
        p_solver_->setX(p_x_);
        p_solver_->setY(p_y_);
        p_solver_->setBeta(&beta_);
        //initialize the solver
        p_solver_->initializeSolver();
      }
      /**initialization of the class with a new beta0
       * @param beta initial start for beta
       * */
      void initializeBeta(STK::VectorX const& beta)
      {
        beta_ = beta;
        p_solver_->initializeSolver();
      }
  };
} // HD


#endif /* LOGISTICLASSO_H_ */
