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

/** @file EnetSolver.cpp
 *  @brief In this file .
 **/

#include "EnetSolver.h"

namespace HD
{
  /*default constructor*/
  EnetSolver::EnetSolver()
  : IPenalizedSolver(), p_solver_(0), Xty_(), b_(), nbActiveVariables_(0), threshold_(), p_penalty_(0)
  {
  }

  /*
   * Constructor
   * @param p_x pointer to the current Data
   * @param beta initial solution
   * @param p_y pointer to the response
   * @param threshold threshold for shrinkage
   * @param p_solver pointer to the solver
   * @param p_penalty pointer to the lasso penalty
   */
  EnetSolver::EnetSolver(STK::ArrayXX const* p_x, STK::VectorX const& beta, STK::VectorX const* p_y,
                           STK::Real const& threshold, STK::CG<EnetMultiplicator,STK::VectorX,InitFunctor>* p_solver, EnetPenalty* p_penalty)
              : IPenalizedSolver(beta,p_x,p_y)
              , p_solver_(p_solver)
              , Xty_()
              , b_()
              , nbActiveVariables_(beta.sizeRows())
              , threshold_(threshold)
              , p_penalty_(p_penalty)
  {
    Xty_ = currentX_.transpose() * *p_y;
    p_solver_->setB(b_);
  }

  /*Initialization of the solver*/
  Real EnetSolver::initializeSolver()
  {
    //check the existence pointers to the data and the response
    if(p_x_ == 0)
      throw STK::invalid_argument(STK::String("p_x_ has not be set"));
    if(p_y_ == 0)
      throw STK::invalid_argument(STK::String("p_y_ has not be set"));

//    currentX.resize(p_x_->sizeRows(),p_x_->sizeCols());
    currentX_ = *p_x_;
    currentBeta_ = beta_;
    Xty_ = currentX_.transpose() * (*p_y_);
    nbActiveVariables_ = p_x_->sizeCols();
    p_solver_->setB(b_);

    currentSet_.resize(nbActiveVariables_);
    for(int i = 1; i <= nbActiveVariables_; i++)
      currentSet_[i] = i;
    return computeLlc();
  }

  /* Computation of the completed loglikelihood*/
  STK::Real EnetSolver::computeLlc()
  {
    STK::Real llc= -( ( (*p_y_ - (currentX_ * currentBeta_) ).square().sum() )/p_penalty_->sigma2() +  p_penalty_->penaltyTerm(currentBeta_))/2;

    return llc;
  }

  /**
   * run the solver of the M step
   * @return
   */
  STK::Real EnetSolver::run(bool burn)
  {
    //compute the b of the linear system Ax=b
    b_.resize(currentSet_.sizeRows());
    for(int i = 1; i <= b_.sizeRows(); i++)
      b_[i]=Xty_[currentSet_[i]];

    //set the b of the linear system Ax=b
    b_ = p_penalty_->multSqrtInvPenalty(b_);

    //run the conjugate gradient
    p_solver_->run();

    //backtransform the solution x to beta
    currentBeta_ = p_penalty_->multSqrtInvPenalty(p_solver_->x());

    //compute llc
    STK::Real llc;
    llc = computeLlc();


    return llc;
  }

  /*run the update of the penalty*/
  void EnetSolver::update(bool toUpdate)
  {
    if(p_penalty_ == 0)
      throw STK::invalid_argument(STK::String("p_penalty_ has not be set"));
    //if burn, check if we can reduce the dimension
    if(toUpdate)
    {
      //thresholding + update currentBeta_ and currentX_;
      updateCurrentBeta();
    }
    else
      beta_ = currentBeta_;

    p_penalty_->update(currentBeta_);
  }

  /*Thresholding of the new estimates : estimated coefficients < threshold_ become 0*/
  void EnetSolver::thresholding()
  {
    for(int i = 1; i <= currentBeta_.sizeRows(); i++)
      if(std::abs(currentBeta_[i]) < threshold_)
      {
        currentBeta_[i] = 0;
        nbActiveVariables_--;
      }
  }


  /* update all the current variables*/
  void EnetSolver::updateCurrentBeta()
  {
    //thresholding of the new estimates
    thresholding();

    //save the new beta
    for(int i = 1; i <= currentBeta_.sizeRows(); i++)
      beta_[currentSet_[i]]=currentBeta_[i];

    if(nbActiveVariables_ != currentBeta_.sizeRows())//if TRUE, at least one variable became 0 in the previous M step
    {
      //update the currentSet
      int pos=1;
      for(int i = 1; i <= currentBeta_.sizeRows(); i++)
      {
        //erase index of the zero estimatee
        if(currentBeta_[i]==0)
        {
          currentSet_.erase(pos);
          pos--;
        }
        pos++;
      }

      //update current beta and current data
      updateCurrentData();
    }
  }

  /*Update the currentBeta_ and currentX_*/
  void EnetSolver::updateCurrentData()
  {
    //resize currentBeta_,  currentX_ with only the active variable
    currentBeta_.resize(nbActiveVariables_);
    currentX_.resize(p_x_->sizeRows(), nbActiveVariables_);

    for(int i = 1; i <= nbActiveVariables_; i++)
    {
      currentBeta_[i]=beta_[currentSet_[i]];
      currentX_.col(i)=p_x_->col(currentSet_[i]);
    }
  }

}


