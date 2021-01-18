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
 * created on: 21 juin 2013
 * Author:   Quentin Grimonprez
 **/

/** @file LassoSolver.h
 *  @brief In this file, Definition of the solver for a lasso penalty .
 **/


#ifndef LASSOSOLVER_H_
#define LASSOSOLVER_H_

#include "IPenalizedSolver.h"
#include "LassoPenalty.h"

namespace HD
{

  /** @ingroup lassoModels
   *  @brief This class inherits from the @c IPenalizedSolver class. It implements the way to solve the Mstep for a lasso penalty
   */
  class LassoSolver : public IPenalizedSolver
  {
    public:
      typedef STK::CG<LassoMultiplicator,STK::VectorX,InitLassoFunctor> CGSolver;

      /** Constructor
       * @param p_x pointer to the current Data
       * @param p_y pointer to the response
       * @param p_beta pointer on the initial/final solution
       * @param threshold threshold for shrinkage
       * @param epsCG tolerance for the Conjugate Gradient
       * @param p_penalty pointer to the lasso penalty
       */
      LassoSolver(LassoPenalty* p_penalty);
      LassoSolver( STK::ArrayXX const* p_x
                 , STK::VectorX const* p_y
                 , STK::VectorX* p_beta
                 , STK::Real const& threshold = 1e-10
                 , STK::Real const& epsCG = 1e-8
                 , LassoPenalty* p_penalty = 0 );
      /**destructor*/
      inline virtual ~LassoSolver() {};
      /**Initialization of the solver (to use when data are modified) */
      STK::Real initializeSolver();
      /** update the solver of the solver (to use when lambda is modified) */
      virtual STK::Real updateSolver();
      /** run the update of the penalty (called during the eStep) */
      void update(bool toUpdate);
      /**Solve the M-step with a conjugate gradient
       * @return the completed loglikelihood
       * */
      STK::Real run(bool toUpdate);
      //getter
      /**@return the pointer to the penalty*/
      inline LassoPenalty* p_penalty() const { return p_penalty_;}
      /**@return a pointer to the CG solver*/
      inline CGSolver const* p_gcsolver() const { return &cgsolver_;}
      /** set the LassoPenalty
       *  @param p_penalty pointer to the penalty
       */
      inline void setPenalty(LassoPenalty* p_penalty) {p_penalty_ = p_penalty;}
      /** @param eps tolerance of the CG */
      inline void setCGEps(STK::Real const& eps) { cgsolver_.setEps(eps); }
      /** Computation of the completed loglikelihood using the current beta */
      STK::Real computeLlc() const;

    protected:
      /** compute an initial value of beta using ols */
      void computeInitialBeta();
      /** update all the current variables*/
      void updateCurrentBeta();
      /** update the system */
      void updateSystem();
      /** update the second member */
      void updateB();

    private:
      ///t(X) * y
      STK::VectorX Xty_;
      ///b from ax=b for CG
      STK::VectorX b_;
      /// x0 for CG
      STK::VectorX x0_;
      ///pointer to the lasso penalty
      LassoPenalty* p_penalty_;
      /// multiplicator for conjugate gradient
      LassoMultiplicator mult_;
      /// conjugate gradient for solver
      CGSolver cgsolver_;
      /// initialization of the  cg solver
      InitLassoFunctor cginit_;
  };
}

#endif /* LASSOSOLVER_H_ */
