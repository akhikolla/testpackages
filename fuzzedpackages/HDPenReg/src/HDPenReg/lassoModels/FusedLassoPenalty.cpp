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
 * created on: 13 sept. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file FusedLassoPenalty.cpp
 *  @brief In this file, implementation of the method of the @c FusedLassoPenalty class .
 **/

#include "FusedLassoPenalty.h"

namespace HD
{
  /* Constructor
   *  @param lambda1 penalization parameter for the l1-norm of the estimates
   *  @param lambda2 penalization parameter for the l1-norm of the difference between successive estimates
   *  @param eps epsilon to add to denominator of fraction to avoid zeros.
   */
  FusedLassoPenalty::FusedLassoPenalty( STK::Real lambda1, STK::Real lambda2, STK::Real eps)
                                      : lambda1_(lambda1)
                                      , lambda2_(lambda2)
                                      , mainDiagonal_()
                                      , offDiagonal_()
                                      , sigma2_(1)
                                      , eps_(eps)
  {}
  /*
   * Copy constructor
   * @param penalty LassoPenalty object to copy
   */
  FusedLassoPenalty::FusedLassoPenalty( FusedLassoPenalty const& penalty)
                                      : IPenalty(penalty)
                                      , lambda1_(penalty.lambda1())
                                      , lambda2_(penalty.lambda2())
                                      , mainDiagonal_(penalty.mainDiagonal())
                                      , offDiagonal_(penalty.offDiagonal())
                                      , sigma2_(penalty.sigma2())
                                      , eps_(penalty.eps())
  {}
  /* clone*/
  FusedLassoPenalty* FusedLassoPenalty::FusedLassoPenalty::clone() const
  { return new FusedLassoPenalty(*this);}

  /* @param beta current estimates
   * @return t(beta) * matrixB * beta
   */
  STK::Real FusedLassoPenalty::penaltyTerm(STK::VectorX const& beta) const
  {
    STK::Real pen = lambda1_ * beta.abs().sum();
    STK::Range r0(beta.begin(), beta.size()-1), r1(beta.begin()+1, beta.size()-1);
    pen += lambda2_ * (beta.sub(r0) - beta.sub(r1)).abs().sum();
    return pen;
  }

  /* update the penalty matrix : matrixB_
   *  @param beta current estimates
   */
  void FusedLassoPenalty::update(STK::VectorX const& beta)
  {
    //resize to the current size
    offDiagonal_.resize(beta.size()-1);
    mainDiagonal_.resize(beta.range());

    if(beta.size() == 1)
    {  mainDiagonal_.front() = lambda1_/(std::abs(beta.front()) + eps_);}
    else
    {
      // first element
      offDiagonal_.front()  = -lambda2_/(std::abs(beta[beta.begin()+1]-beta.front()) + eps_);
      mainDiagonal_.front() =  lambda1_/(std::abs(beta.front()) + eps_) - offDiagonal_.front();
      // intermediate elements
      for(int i = beta.begin()+1; i < beta.lastIdx(); i++)
      {
        offDiagonal_[i] = -lambda2_/(std::abs(beta[i+1]-beta[i]) + eps_);
        mainDiagonal_[i] = lambda1_/(std::abs(beta[i]) + eps_) - offDiagonal_[i] - offDiagonal_[i-1];
      }
      // last element
      mainDiagonal_[beta.lastIdx()] = lambda1_/(std::abs(beta[beta.lastIdx()]) + eps_) - offDiagonal_[beta.lastIdx()-1];
    }
  }

  /* update the penalty matrix : matrixB_
   *  @param beta current estimates
   *  @param segment segment repartition from FusedLassoSolver
   */
  void FusedLassoPenalty::update(STK::VectorX const& beta, Segment const& segment)
  {
    //resize to the current size
    offDiagonal_.resize(beta.size()-1);
    mainDiagonal_.resize(beta.range());

    //we have to multiply lambda1 by the number of element of the segment to keep the tight structure when regrouping variable
    if(beta.size() == 1)
      mainDiagonal_.front() = segment.front().size() * lambda1_/(std::abs(beta.front()) + eps_);
    else
    {
      // first element
      offDiagonal_.front()  = -lambda2_/(std::abs(beta[beta.begin()+1]-beta.front()) + eps_);
      mainDiagonal_.front() = segment.front().size() * lambda1_/(std::abs(beta.front()) + eps_) - offDiagonal_.front();
      // intermediate elements
      for(int i = beta.begin()+1; i < beta.lastIdx(); i++)
      {
        offDiagonal_[i] = -lambda2_/(std::abs(beta[i+1]-beta[i]) + eps_);
        mainDiagonal_[i] = segment[i].size() * lambda1_/(std::abs(beta[i]) + eps_) - offDiagonal_[i] - offDiagonal_[i-1];
      }
      // last element
      mainDiagonal_.back() = segment.back().size() *  lambda1_/(std::abs(beta.back()) + eps_) - offDiagonal_[beta.lastIdx()-1];
    }
  }

}


