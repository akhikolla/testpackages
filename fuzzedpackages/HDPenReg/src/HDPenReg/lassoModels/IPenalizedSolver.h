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
 * created on: 14 juin 2013
 * Author:   Quentin Grimonprez
 **/

/** @file IPenalizedSolver.h
 *  @brief In this file definition of the interface class @c IPenalizedSolver.
 **/


#ifndef IPENALIZEDSOLVER_H_
#define IPENALIZEDSOLVER_H_

#include "IPenalty.h"

namespace HD
{
  /**functor for Initial CG  */
  struct InitLassoMultiplicator
  {
      /** functor for the CG. A=X'*X
       * @param x vector of length sizeCols(A)
       * @return A*x
       */
      STK::VectorX operator()(STK::VectorX const& x) const
      {
        STK::VectorX a = p_x_->transpose() * (*p_x_  * x) + lambda_ * x;
        return a;
      }
      /** Constructor of the functor
       * @param p_x constant pointer on the data
       */
      InitLassoMultiplicator( STK::ArrayXX const* p_x, STK::Real const& lambda): p_x_(p_x), lambda_(1.) {}
      ///pointer to the current data
      STK::ArrayXX const* p_x_;
      STK::Real lambda_;
  };

  /** Functor for initialization in the conjugate gradient */
  struct InitLassoFunctor
  {
      /** @return pointer on the initial value */
      inline STK::VectorX operator()() const  { return *p_x0_;}
      /** Constructor
       *  @param p_x pointer on the initial value
       **/
      inline InitLassoFunctor() : p_x0_(0) {};
      ///pointer on the initial value
      STK::VectorX const* p_x0_;
  };

  /** Functor for initialization in the conjugate gradient */
  struct InitFunctor
  {
      /** @return pointer on the initial value */
      inline STK::VectorX operator()() const  { return *p_x_*lambda_;}
      /** Constructor
       *  @param p_x pointer on the initial value
       **/
      inline InitFunctor(STK::VectorX const* p_x = 0) : p_x_(p_x), lambda_(1.) {};
      ///pointer on the initial value
      STK::VectorX const* p_x_;
      STK::Real lambda_;
  };

  /** Functor for initialization in the conjugate gradient */
  struct InitIdFunctor
  {
      /** @return pointer on the initial value */
      inline STK::VectorX const& operator()() const  { return *p_x_;}
      /** Constructor
       *  @param p_x pointer on the initial value
       **/
      inline InitIdFunctor(STK::VectorX const* p_x = 0) : p_x_(p_x) {};
      ///pointer on the initial value
      STK::VectorX const* p_x_;
  };

  /** @ingroup lassoModels
   *  @brief The class IPenalizedSolver is an interface for the solver
   *  of the @c PenalizedModels M-step
   */
  class IPenalizedSolver
  {
    public:
      /**default constructor*/
      IPenalizedSolver(): currentX_()
                        , currentBeta_()
                        , currentSet_()
                        , p_beta_(0)
                        , p_x_(0)
                        , p_y_(0)
                        , threshold_(1e-10)
      {}

      /** Constructor
       *  @param beta value for initializing beta
       *  @param p_x pointer to the data
       *  @param p_y pointer to the response
       */
      IPenalizedSolver( STK::VectorX* p_beta
                      , STK::ArrayXX const* p_x
                      , STK::VectorX const* p_y
                      , STK::Real threshold = 1e-10
                      )
                      : currentX_(*p_x)
                      , currentBeta_(*p_beta)
                      , currentSet_(p_x->sizeCols())
                      , p_beta_(p_beta)
                      , p_x_(p_x)
                      , p_y_(p_y)
                      , threshold_(threshold)
       {
         for(int i = currentSet_.begin(); i < currentSet_.end(); i++)
           currentSet_[i] = i;
       };
       /**destructor*/
       inline virtual ~IPenalizedSolver() {};

       /**run the update of the penalty (Estep)*/
       virtual void update(bool toUpdate) = 0;
       /** run the solver (Mstep)
        *  @return the new estimated value of beta
        */
       virtual STK::Real run(bool burn) = 0;
       /** initialize the initial beta0 (do nothing by default)*/
       virtual STK::Real updateSolver() = 0;

       /** initialize all the containers of the class */
       virtual STK::Real initializeSolver() = 0;
       /**@return the current threshold */
       inline STK::Real const& threshold() const { return threshold_; }

       //setter
       inline void setThreshold(STK::Real const& threshold) { threshold_ = threshold; }
       /** set the pointer to the current data
        *  (data reduced to covariates from current set)*/
       inline void setX(STK::ArrayXX const* p_x) { p_x_ = p_x;};
       /**set the pointer to the response*/
       inline void setY(STK::VectorX const* p_y) { p_y_ = p_y;};
       /**set the pointer to the current beta*/
       inline void setBeta(STK::VectorX* p_beta)
       {
         p_beta_ = p_beta;
         currentBeta_ = *p_beta;
       };

       //getter
       /**@return beta_*/
       inline STK::VectorX* p_beta() const { return p_beta_;};
       /**@return currentX_*/
       inline STK::ArrayXX const* p_currentX() const {return &currentX_;};
       /**@return currentX_*/
       inline STK::ArrayXX const& currentX() const {return currentX_;};
       /**@return currentBeta_*/
       inline STK::VectorX const& currentBeta() const {return currentBeta_;};
       ///@return currentSet_
       inline STK::VectorXi const& currentSet() const {return currentSet_;};
       /**@return currentBeta_*/
       inline STK::VectorX const* p_currentBeta() const {return &currentBeta_;};

       /** disrupts current value of beta */
       void disruptsBeta()
       {
         for (int i=p_beta_->begin(); i< p_beta_->end(); ++i)
         {
           // disrupts only small values
           if (std::abs(p_beta_->elt(i))< threshold_)
           {
             p_beta_->elt(i) += STK::sign(p_beta_->elt(i), 10. * threshold_);
           }
         }
       }

    protected:
       /** compute the loglikelihood
        * @return the loglikelihood of the current step
        */
       virtual STK::Real computeLlc() const = 0;

    protected:
       /// current data
       STK::ArrayXX currentX_;
       /// Current beta
       STK::VectorX currentBeta_;
       /// current set of variable
       STK::VectorXi currentSet_;
       /// pointer on beta_
       STK::VectorX* p_beta_;
       /// pointer on the data
       STK::ArrayXX const* p_x_;
       /// pointer on the response
       STK::VectorX const* p_y_;
       ///threshold under we consider a beta equal to 0
       STK::Real threshold_;
   };
}

#endif /* IPENALIZEDSOLVER_H_ */
