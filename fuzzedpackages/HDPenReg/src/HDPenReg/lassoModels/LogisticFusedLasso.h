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

/** @file LogisticFusedLasso.h
 *  @brief In this file .
 **/


#ifndef LOGISTICFUSEDLASSO_H_
#define LOGISTICFUSEDLASSO_H_

#include "PenalizedModels.h"
#include "FusedLassoPenalty.h"
#include "LogisticFusedLassoSolver.h"

namespace HD
{
  class LogisticFusedLasso;

  template<>
  struct ModelTraits<LogisticFusedLasso>
  {
    typedef LogisticFusedLassoSolver Solver;
    typedef FusedLassoPenalty Penalty;
    typedef FusedLassoMultiplicator Multiplicator;
    typedef STK::CG<FusedLassoMultiplicator,STK::VectorX, InitIdFunctor> CG;
  };

  /**
   * Class FusedLasso derived from @c PenalizedModels.
   * This class constructs a FusedLasso Model to be solved by an @c EM algorithm.
   */
  class LogisticFusedLasso : public PenalizedModels<LogisticFusedLasso>
  {
    public:
      typedef STK::CG<FusedLassoMultiplicator,STK::VectorX,InitFunctor> CGSolver;
      /** constructor
       *  @param p_x pointer to the data
       *  @param p_y pointer to the response
       *  @param lambda1 value of parameter associated to the l1 penalty
       *  @param lambda2 value of parameter associated to the l1 penalty of successive coefficients
       *  @param threshold threshold for setting coefficient to 0
       *  @param epsCG epsilon for CG convergence
       */
      LogisticFusedLasso( STK::ArrayXX const* p_x, STK::VectorX const* p_y
                        , STK::Real lambda1, STK::Real lambda2
                        , STK::Real threshold, STK::Real epsCG)
                        : PenalizedModels<LogisticFusedLasso>(p_x, p_y)
      {
#ifdef HD_DEBUG
  std::cout << "Creating LogisticFusedLasso::LogisticFusedLasso with"
            << "lambda1 =" << lambda1 << "\n"
            << "lambda2 =" << lambda2 << "\n"
            << "threshold =" << threshold << "\n"
            << "epsCGC =" << epsCG << std::endl;
#endif
  // creation lasso penalty
  p_penalty_ = new FusedLassoPenalty(lambda1, lambda2, 1e-08);
  //create solver for lasso and add the penalty
  p_solver_ = new LogisticFusedLassoSolver(p_x_, p_y_, &beta_, threshold, epsCG, p_penalty_);
  // creation lasso penalty
#ifdef HD_DEBUG
  std::cout << "LogisticFusedLasso Initialized. Likelihood =" << p_solver_->computeLlc() << std::endl;
#endif
      }

      /** destructor*/
      inline ~LogisticFusedLasso() {}
      /** set the lasso regularization parameter
       * @param lambda1
       */
      inline void setLambda1(STK::Real const& lambda1) {p_penalty_->setLambda1(lambda1);}
      /** set the fusion regularization parameter
       * @param lambda2
       */
      inline void setLambda2(STK::Real const& lambda2) {p_penalty_->setLambda2(lambda2);}
      /** set the threshold for segments
       * @param threshold
       */
      inline void setThreshold(STK::Real const& threshold) {p_solver_->setThreshold(threshold);}
      /** set the epsilon for avoid 0 to denominator
       * @param epsilon
       */
      inline void setEps(STK::Real const& eps) {p_penalty_->setEps(eps);}
      /** set the epsilon for the CG
       *  @param epsilon epsilon for the convergence of CG
       */
      inline void setCGEps( STK::Real eps) { p_solver_->setCGEps(eps);}

      /**initialize the containers of all subclasses*/
      void initializeModel()
      {
        //set the parameters of the FusedLassoSolver
        p_solver_->setX(p_x_);
        p_solver_->setY(p_y_);
        p_solver_->setBeta(&beta_);
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
}


#endif /* LOGISTICFUSEDLASSO_H_ */
