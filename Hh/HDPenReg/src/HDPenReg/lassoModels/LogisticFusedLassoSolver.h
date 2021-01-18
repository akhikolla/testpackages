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

/** @file LogisticFusedLassoSolver.h
 *  @brief In this file .
 **/


#ifndef LOGISTICFUSEDLASSOSOLVER_H_
#define LOGISTICFUSEDLASSOSOLVER_H_


#include "IPenalizedSolver.h"
#include "FusedLassoPenalty.h"

namespace HD
{


  /** @ingroup lassoModels
   *  @brief This class inherits from the @c IPenalizedSolver class. It implements the way to solve the Mstep for a Fused lasso penalty
   */
  class LogisticFusedLassoSolver : public IPenalizedSolver
  {
    public:
      typedef STK::CG<FusedLassoMultiplicator, STK::VectorX,InitIdFunctor> CGSolver;
      typedef STK::Array2DVector<STK::Range> Segment;

      /**
       * Constructor
       * @param p_x pointer to the full data
       * @param beta initial solution of the problem
       * @param p_y pointer to the response associated with the data
       * @param p_cgsolver pointer to the solver
       * @param p_penalty pointer to the Fused lasso penalty
       * @param eps tolerance for zero
       */
      LogisticFusedLassoSolver( STK::ArrayXX const* p_x, STK::VectorX const* p_y, STK::VectorX* beta
                              , STK::Real const& threshold, STK::Real const& epsCG
                              , FusedLassoPenalty* p_penalty);

      /**destructor*/
      inline ~LogisticFusedLassoSolver() {};

      //getter
      /**@return A constant pointer on z*/
      inline STK::VectorX const*  p_z() const { return &z_;}
      /**@return  z*/
      inline STK::VectorX const&  z() const { return z_;}
      /** @return segment_ */
      inline Segment segment() const {return segment_;}
      /**@return pointer to the penalty*/
      inline FusedLassoPenalty* p_penalty() {return p_penalty_;}
      /**@return pointer to the solver*/
      inline CGSolver* p_solver() { return &cgsolver_;}

      //setter
      /** set the pointer to the FusedLassopenalty
       * @param p_penalty pointer to the FusedLassopenalty
       * */
      inline void setPenalty(FusedLassoPenalty* p_penalty) {p_penalty_ = p_penalty;}
      /** set the threshold
       *  @param threshold threshold for shrinkage to 0
       */
      inline void setThreshold(STK::Real threshold) { threshold_ = threshold;}
      /** set the tolerance for the difference of estimates
       * @param eps new tolerance
       */
      inline void setEps(STK::Real const& eps) {eps_ = eps;}
      /** @param eps tolerance of the CG */
      inline void setCGEps(STK::Real const& eps) { cgsolver_.setEps(eps); }

      /** Solve the M-step with a conjugate gradient
       *  @return the completed loglikelihood
       */
      STK::Real run(bool toUpdate);
      /** run the update of the penalty*/
      void update(bool toUpdate);

      /**Initialization of the solver*/
      STK::Real initializeSolver();
      /** initialize beta0 */
      virtual STK::Real updateSolver();
      /**Computation of the completed log-likelihood
       * @return the current completed loglikelihood
       */
      STK::Real computeLlc() const;

    protected:
      /** compute an initial value of beta, and set it to currentBeta_ */
      void computeInitialBeta();
      /** update currentSet_, segment_ and currentBeta_
       *  @return a boolean, if true there is a changement in the currentSet
       */
      bool updateCurrentBeta();
      /** update the system */
      void updateSystem();
      /** update Z values using Expectation */
      void updateZ();
      /**update currentX_*/
      void updateCurrentData();
      /** update beta_*/
      void updateBeta();

    private:
      ///estimated response z
      STK::VectorX z_;
      ///b from ax=b for CG
      STK::VectorX b_;
      /// vector containing the segment corresponding to every index of beta
      Segment segment_;
      ///size of currentBeta_
      int nbActiveVariables_;
      ///eps_ tolerance for the difference of estimates
      STK::Real eps_;
      /// gaussian law
      STK::Law::Normal normal_;

      ///pointer to the lasso penalty
      FusedLassoPenalty* p_penalty_;

      /// multiplicator for conjugate gradient
      FusedLassoMultiplicator mult_;
      /// conjugate gradient for solver
      STK::CG<FusedLassoMultiplicator,STK::VectorX, InitIdFunctor> cgsolver_;
      /// initial functor for CG
      InitIdFunctor cginit_;
  };
}


#endif /* LOGISTICFUSEDLASSOSOLVER_H_ */
