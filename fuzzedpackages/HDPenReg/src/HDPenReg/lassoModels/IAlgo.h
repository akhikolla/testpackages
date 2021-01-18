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
 * created on: 23 mai 2013
 * Author:   Quentin Grimonprez
 **/

/** @file IAlgo.h
 *  @brief In this file, definition of the interface class @c IAlgo.
 **/


#ifndef IALGO_H_
#define IALGO_H_

#include "PenalizedModels.h"

namespace HD
{
  /** Interface class for algorithms.
   *  Inherited classes has to implemented a run function.
   */
  class IAlgo
  {
    public:
      /**default constructor*/
      IAlgo() : maxStep_(0), eps_(0.) {};
      /** constructor*/
      IAlgo( int maxStep, STK::Real const& eps )
           : maxStep_(maxStep), eps_(eps) {};

      /**virtual destructor*/
      virtual ~IAlgo() {};
//      /**
//       * run an algorithm to solve the penalizedModels
//       * @param model a pointer to the penalizedModels object to solve.
//       */
//      virtual bool run(PenalizedModelsBase* model) = 0;

      /** set the threshold for convergence of the complete loglikelihood
       * @param eps new value for the threshold
       */
      inline void setEps(STK::Real eps) {eps_=eps;}
      /** set the maximum number of step
       * @param maxStep new maximal number of step
       */
      inline void setMaxStep(int maxStep) {maxStep_=maxStep;}

      /** @return the number of step*/
      inline STK::Real eps() const {return eps_;}
      /** @return the maximum number of steps*/
      inline int maxStep() const {return maxStep_;}

    protected:
      ///maximum number of step
      int maxStep_;
      ///threshold for convergence of the complete loglikelihood
      STK::Real eps_;
  };
}

#endif /* IALGO_H_ */
