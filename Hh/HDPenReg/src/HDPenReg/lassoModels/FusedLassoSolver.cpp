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
 * created on: 18 sept. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file FusedLassoSolver.cpp
 *  @brief In this file, implementation of the methods of the @c FusedLassoSolver class. .
 **/

#include "FusedLassoSolver.h"

namespace HD
{
  /*
   * Constructor
   * @param p_x pointer to the full data
   * @param beta initial solution of the problem
   * @param p_y pointer to the response associated with the data
   * @param burn burn-in period before regrouping variables in segments
   * @param p_solver pointer to the solver
   * @param p_penalty pointer to the Fused lasso penalty
   * @param eps tolerance for zero
   */
FusedLassoSolver::FusedLassoSolver( STK::ArrayXX const* p_x, STK::VectorX const* p_y, STK::VectorX* p_beta
                                  , STK::Real const& threshold, STK::Real const& epsCG
                                  , FusedLassoPenalty* p_penalty)
                                  : IPenalizedSolver(p_beta, p_x, p_y, threshold)
                                  , currentXty_()
                                  , segment_(p_x_->cols())
                                  , nbActiveVariables_(p_x_->sizeCols())
                                  , eps_(threshold)
                                  , p_penalty_(p_penalty)
  {
#ifdef HD_DEBUG
    std::cout << "Entering FusedLassoSolver::FusedLassoSolver()\n";
    std::cout << "threshold_ =" << threshold_ << "\n";
    std::cout << "eps_ =" << eps_ << "\n";
    std::cout << "epsCG =" << epsCG << "\n";
#endif
    // currentX_  is initialized to *p_x in IPenalizedSolver
    // currentSet_ is initialized in IPenalizedSolver to currentSet_[i} = i
    // currentXty_ and currentBeta_ are initialized in computeInitialBeta()
    computeInitialBeta();
    // initialize size 1 segments
    for(int i = segment_.begin(); i < segment_.end(); i++)
    { segment_[i] = STK::Range(i,1);}
    // initialize CG
    cgsolver_.setMultFunctor(mult_);
    cgsolver_.setEps(epsCG);
    cgsolver_.setB(currentXty_);
    //intialize mult functor for CG
    mult_.p_x_            = p_currentX();
    mult_.p_mainDiagonal_ = p_penalty_->p_mainDiagonal();
    mult_.p_offDiagonal_  = p_penalty_->p_offDiagonal();
    mult_.p_sigma2_       = p_penalty_->p_sigma2();
    // initialize init for CG
    init_.p_x_ = p_currentBeta();
    cgsolver_.setInitFunctor(&init_);
    p_penalty_->update(currentBeta_);
  }
  /*initialize the container of the class*/
  STK::Real FusedLassoSolver::initializeSolver()
  {
#ifdef HD_DEBUG
    if(p_x_ == 0)
      throw STK::invalid_argument(STK::String("p_x_ has not be set"));
    if(p_y_ == 0)
      throw STK::invalid_argument(STK::String("p_y_ has not be set"));
#endif
    currentX_ = *p_x_;
    currentBeta_ = *p_beta_;
    nbActiveVariables_ = p_x_->sizeCols();
    currentXty_ = currentX_.transpose() * (*p_y_);
    // init sets and segments
    currentSet_.resize(nbActiveVariables_);
    segment_.resize(nbActiveVariables_);
    for(int i = segment_.begin(); i < segment_.end(); i++)
    {
      segment_[i] = STK::Range(i,1);
      currentSet_[i] =  i;
    }
    p_penalty_->update(currentBeta_);
    return updateSolver();
  }
  /*Initialization of the solver*/
  STK::Real FusedLassoSolver::updateSolver()
  {
#ifdef HD_DEBUG
    //check the existence pointers to the data and the response
    if(p_x_ == 0)
      throw STK::invalid_argument(STK::String("p_x_ has not be set"));
    if(p_y_ == 0)
      throw STK::invalid_argument(STK::String("p_y_ has not be set"));
#endif
    // reinitialize all containers
    currentX_ = *p_x_;
    currentBeta_ = *p_beta_;
    nbActiveVariables_ = p_x_->sizeCols();
    currentXty_ = currentX_.transpose() * (*p_y_);
    // init sets and segments
    currentSet_.resize(nbActiveVariables_);
    segment_.resize(nbActiveVariables_);
    for(int i = segment_.begin(); i < segment_.end(); i++)
    {
      segment_[i] = STK::Range(i,1);
      currentSet_[i] =  i;
    }
    p_penalty_->update(currentBeta_);
#ifdef HD_DEBUG
      std::cout << "FusedLassoSolver::updateSolver done. llc =" << computeLlc() << std::endl;
#endif
    return computeLlc();
  }


  /*Computation of the completed log-likelihood
   * @return the current completed loglikelihood
   */
  STK::Real FusedLassoSolver::computeLlc() const
  {
    return - ( ( (*p_y_ - (currentX_ * currentBeta_) ).square().sum() )/p_penalty_->sigma2()
               +  p_penalty_->penaltyTerm(currentBeta_))/2;
  }
  /*update the penalty (E step)*/
  void FusedLassoSolver::update(bool toUpdate)
  {
#ifdef HD_DEBUG
    if(p_penalty_ == 0)
      throw STK::invalid_argument(STK::String("p_penalty_ has not be set"));
#endif
    //reduction of the data with segments
    if (toUpdate)
    {
      //update segment_ and currentSet_ and currentBeta_
      bool changement = updateCurrentBeta();
      //update the full beta_
      updateBeta();
      //if true, we have to change the currentX matrix
      if(changement)
      {
        updateCurrentData();
        //update currentXty_ because currentX_ change
        currentXty_ = currentX_.transpose() * (*p_y_);
      }
    }
    else
    {
      if(p_beta_->range()!=currentBeta_.range())
      { updateBeta();}
      else
      { *p_beta_=currentBeta_;}
    }
    p_penalty_->update(currentBeta_,segment_);
  }
  /* Solve the M-step with a conjugate gradient
   * @return new estimate of beta
   */
  STK::Real FusedLassoSolver::run(bool toUpdate)
  {
    //run the conjugate gradient
    int cgiter = cgsolver_.run();
#ifdef HD_DEBUG
      std::cout << "In FusedLassoSolver::run. cgsolver_ run in " << cgiter << " iterations.\n";
#endif
    //get the solution
    currentBeta_ = cgsolver_.x();
    return computeLlc();
  }

  /* compute an initial value of beta using ols */
  void FusedLassoSolver::computeInitialBeta()
  {
    // compute initial value using OLS
    currentXty_ = currentX_.transpose() * *p_y_;
    InitLassoMultiplicator mult(p_x_, p_penalty_->lambda1());
    STK::CG<InitLassoMultiplicator,STK::VectorX,InitFunctor> cginitial(mult, currentXty_, 0, 1e-05);
    cginitial.run();
    p_beta_->move(cginitial.x());
    disruptsBeta();
    currentBeta_ = *p_beta_;
  }
  /*
   * update currentSet_, segment_ and currentBeta_
   * @return a boolean, if true there is a changement in the currentSet
   */
  bool FusedLassoSolver::updateCurrentBeta()
  {
#ifdef HD_DEBUG
    std::cout << "Entering FusedLassoSolver::updateCurrentBeta()\n";
    std::cout << "p_beta_->range() =" << p_beta_->range() << "\n";
    std::cout << "currentSet_.range() =" << currentSet_.range() << "\n";
    std::cout << "currentBeta_.range() =" << currentBeta_.range() << "\n";
    std::cout << "segment_.range() =" << segment_.range() << "\n";
#endif
    bool changement = false;
    // start fusion of the segment
    for(int i = currentBeta_.lastIdx(); i > currentBeta_.begin(); i--)
    {
      // TRUE: fusion of 2 segments
      if( std::abs(currentBeta_[i]-currentBeta_[i-1]) <= eps_ )
      {
        //increase the last index of the segment i-1
        segment_[i-1] = segment_[i-1].incLast(segment_[i].size());
        //delete the segment i
        //erase the i-th beta (same as (i-1)-th)
        segment_.erase(i);
        currentBeta_.erase(i);
        //the segment i merged with segment i-1, segment i+1, i+2,... become segment i, i+1, ...
        for( int j = i; j < currentSet_.end(); j++) currentSet_[j]--;
        // a variable is removed from the active variable
        nbActiveVariables_--;
        changement = true;
      }
    }
    // check small values
    for(int i = currentBeta_.begin(); i < currentBeta_.end(); i++)
    { if( std::abs(currentBeta_[i]) < eps_ ) currentBeta_[i] = 0;}
    return changement;
  }

  /*
   * update currentX_
   */
  void FusedLassoSolver::updateCurrentData()
  {
    //resize currentX_ with only the active variable
    currentX_.resize(p_x_->rows(), nbActiveVariables_);
    currentX_.zeros();
    for(int i = currentX_.beginCols(); i < currentX_.endCols(); i++)
    {
      //the covariate associates to an index is the sum of all covariates of the segment
      for(int j = segment_[i].begin(); j< segment_[i].end(); j++)
        currentX_.col(i) += p_x_->col(j);
    }
  }

  /*
   * update beta_
   */
  void FusedLassoSolver::updateBeta()
  {
    //complete beta with the updated values
    for(int i = currentBeta_.begin(); i < currentBeta_.end(); i++)
    {
      for(int j = segment_[i].begin(); j < segment_[i].end(); j++)
        p_beta_->elt(j) = currentBeta_[i];
    }
  }

}

