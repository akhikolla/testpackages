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
 * created on: 3 avr. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file Fusion.h
 *  @brief In this file, we define the @c Fusion class.
 **/


#ifndef FUSION_H_
#define FUSION_H_

namespace HD
{
/**
 * This class transform the input data in order to use the lars algorithm for the fusion problem.
 *

 *
 * Let \f$ X\f$ a matrix of size \f$ n\times p\f$, \f$ y\f$ a vector of length \f$ n\f$ and \f$\lambda\f$ a positive real.
 *
 * The fusion problem is to find \f$ \hat{\beta}\f$ such that
 *
 * \f$ \hat{\beta}=argmin_{\beta} \|\mathbf{y}-\mathbf{X}\mathbf{\beta}\|_2^2 + \lambda\sum\limits_{i=1}^{p-1} |\beta_{i+1}-\beta_i|\f$
 *
 *
 * This problem can be solve with the LARS algorithm by transforming it into a lasso problem.
 *
 * Let \f$ L^-1\f$ a lower triangular matrix of 1. Transforming \f$ X \f$ to \f$ Z=XL^{-1}\f$ and \f$ \theta=L\beta\f$ and
 * you have a lasso problem in \f$\theta\f$.
 *
 *
 */
  class Fusion
  {
    public:
      //constructors
      /**
       * constructor
       * @param X matrix of size n*p, a row contains the values of each covariate for an individual.
       * @param y vector of length n containing the response
       * @param intercept if true, there is an intercept in the model
       */
      Fusion(STK::CArrayXX const& X, STK::CVectorX const& y, bool intercept = true);

      /**
       *
       * @param X matrix of size n*p, a row contains the values of each covariate for an individual.
       * @param y vector of length n containing the response
       * @param maxSteps number of maximum step to do
       * @param intercept if true, there is an intercept in the model
       * @param eps epsilon (for 0)
       */
      Fusion(STK::CArrayXX const& X, STK::CVectorX const& y, int maxSteps, bool intercept = true, STK::Real eps =STK::Arithmetic<STK::Real>::epsilon());


      //getters
      /**@return path of the coefficients*/
      inline Path path() const {return path_;}

      /**@return Number of step of the algorithm*/
      inline int step() const {return step_;}

      /**
       * @param i step
       * @return the Pathstate object : the state of the path at the step i
       */
      inline PathState const& path(int i) const {return path_.states(i);}

      /**
       * @param i index of the step
       * @param j index of the coefficients
       * @return the value of the j-th coefficient at the step i
       */
      inline STK::Real coefficient(int i,int j) const {return path_.varCoeff(i,j);}

      /**
       * @param i index of the step
       * @param j index of the coefficients
       * @return the value of the j-th coefficient at the step i
       */
      inline int varIdx(int i,int j) const {return path_.varIdx(i,j);}

      /**
       * @param i index of the step
       * @return the value of l1norm at the i-th step
       */
      inline STK::Real l1norm(int i) const {return path_.l1norm(i);}

      /**@return the value of lambda */
      inline std::vector<STK::Real> lambda() const {return path_.lambda();}

      /** @return the historic of add and drop variable*/
      inline std::vector< std::pair<std::vector<int>,std::vector<int> > > evolution() const {return path_.evolution();}

      /** @return the intercept of the solution*/
      inline STK::Real mu() const {return mu_;}

      /**@return muX_*/
      inline STK::CVectorX muX() const {return muX_;}

      /**@return muX_[i]*/
      inline STK::Real muX(int i) const {return muX_[i];}

      /** @return the ignored variable*/
      inline STK::CArrayVector<bool> toIgnore() const {return toIgnore_;}

      /** @return msg_error_*/
      inline std::string msg_error() const {return msg_error_;}

      /**
       * run the lars algorithm for solving the fusion problem on Z=X*L^-1 (L^-1 = lower triangular matrix of 1)
       */
      void run();

    protected:
      /**
       * change X in Z=X*L^-1 (L^-1 = lower triangular matrix of 1)
       */
      void computeZ();

    private:
      ///matrix of size n*p, a col = a covariate
      STK::CArrayXX X_;
      ///vector size n, response
      STK::CVectorX y_;
      ///maximum number of steps for the lars algorithm
      int maxSteps_;
      ///numerical zero
      STK::Real eps_;
      ///mean
      STK::Real mu_;
      ///meanX
      STK::CVectorX muX_;
      ///number of step done by the lars algorithm
      int step_;
      ///solution path of the lars
      Path path_;
      ///ignored variables (due to correlation)
      STK::CArrayVector<bool> toIgnore_;
      /// if true, there is an intercept in the model
      bool intercept_;
      ///last error message
      std::string msg_error_;

  };
}


#endif /* FUSION_H_ */
