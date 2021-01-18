/*--------------------------------------------------------------------*/
/*     Copyright (C) 2013-2013  Serge Iovleff, Quentin Grimonprez

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
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
 * created on: 2 oct. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file Lasso.h
 *  @brief In this file, we define the @c Lasso class .
 **/


#ifndef LASSO_H_
#define LASSO_H_

#include "PenalizedModels.h"
#include "LassoPenalty.h"
#include "LassoSolver.h"

namespace HD
{
  class Lasso;

  template<>
  struct ModelTraits<Lasso>
  {
    typedef LassoSolver Solver;
    typedef LassoPenalty Penalty;
    typedef LassoMultiplicator Multiplicator;
    typedef STK::CG<LassoMultiplicator,STK::VectorX, InitFunctor> CG;
  };

  /**
   * Class Lasso derived from @c PenalizedModels.
   * This class constructs a Lasso Model to be solved by an @c EM algorithm.
   */
  class Lasso : public PenalizedModels<Lasso>
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
      //add
      Lasso() : PenalizedModels<Lasso>()
      {
          // creation lasso penalty
          p_penalty_ = new LassoPenalty(0.1);
          //create solver for lasso and add the penalty
          p_solver_ = new LassoSolver(p_penalty_);
      }
      Lasso( STK::ArrayXX const* p_x
           , STK::VectorX const* p_y
           , STK::Real lambda
           , STK::Real threshold, STK::Real epsCG
           )
           : PenalizedModels<Lasso>(p_x, p_y)
      {
#ifdef HD_DEBUG
        std::cout << "Creating Lasso. Lambda =" << lambda << std::endl;
#endif
        // creation lasso penalty
        p_penalty_ = new LassoPenalty(lambda);
        //create solver for lasso and add the penalty
        p_solver_ = new LassoSolver(p_x_, p_y_, &beta_, threshold, epsCG, p_penalty_);
#ifdef HD_DEBUG
        std::cout << "Lasso Initialized. Likelihood =" << p_solver_->computeLlc() << std::endl;
#endif
      }

      /** destructor*/
      virtual ~Lasso() {}
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
      /** initialization of the class using the current value of beta */
      void initializeBeta()
      {
        p_solver_->disruptsBeta();
        p_solver_->initializeSolver();
#ifdef HD_DEBUG
        std::cout << "Lasso::initializeBeta done. Likelihood =" << p_solver_->computeLlc() << std::endl;
#endif
      }
      /**initialize the containers of all subclasses*/
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
}

#endif /* LASSO_H_ */
