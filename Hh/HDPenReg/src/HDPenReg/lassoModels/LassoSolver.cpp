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

/** @file LassoSolver.cpp
 *  @brief In this file, implementation of the methods of the @c LassoSolver class .
 **/

#include "LassoSolver.h"

namespace HD
{
  /*
   * Constructor
   * @param p_x pointer to the current Data
   * @param beta initial solution
   * @param p_y pointer to the response
   * @param threshold threshold for shrinkage
   * @param p_solver pointer to the solver
   * @param p_penalty pointer to the lasso penalty
   */
/** Default Constructor       */
	LassoSolver::LassoSolver(LassoPenalty* p_penalty): IPenalizedSolver()
                            , Xty_()
                            , b_(), x0_()
                            , p_penalty_(p_penalty)
                            , mult_(), cgsolver_(), cginit_()
{}

  LassoSolver::LassoSolver( STK::ArrayXX const* p_x, STK::VectorX const* p_y, STK::VectorX* p_beta
                          , STK::Real const& threshold, STK::Real const& epsCG
                          , LassoPenalty* p_penalty)
                          : IPenalizedSolver(p_beta, p_x, p_y, threshold)
                          , Xty_()
                          , b_(), x0_()
                          , p_penalty_(p_penalty)
                          , mult_(), cgsolver_(), cginit_()
  {
    // currentX_  is initialized to *p_x in IPenalizedSolver
    // currentSet_ is initialized in IPenalizedSolver to currentSet_[i} = i
    // Xty_ and currentBeta_ are initialized in computeInitialBeta
    computeInitialBeta();
    // initialize CG
    cgsolver_.setMultFunctor(mult_);
    cgsolver_.setEps(epsCG);
    cgsolver_.setB(b_);
    cginit_.p_x0_    = &x0_;
    cgsolver_.setInitFunctor(&cginit_);
    //intialize mult functor for CG
    mult_.p_x_              = &currentX_;
    mult_.p_sigma2_         = p_penalty_->p_sigma2();
    mult_.p_sqrtInvPenalty_ = p_penalty_->p_sqrtInvPenalty();
    // initialize penalty term
    p_penalty_->update(currentBeta_);
    x0_ = p_penalty_->sqrtInvPenalty() * Xty_;
  }
  /*Initialization of the solver*/
  STK::Real LassoSolver::initializeSolver()
  {
#ifdef HD_DEBUG
    //check the existence pointers to the data and the response
    if(p_x_ == 0)
      throw STK::invalid_argument(STK::String("p_x_ has not be set"));
    if(p_y_ == 0)
      throw STK::invalid_argument(STK::String("p_y_ has not be set"));
#endif
#ifdef HD_DEBUG
      std::cout << "LassoSolver::initializeSolver done. llc =" << computeLlc() << std::endl;
#endif
      // reinitialize all containers
      currentX_    = *p_x_;
      currentBeta_ = *p_beta_;
      Xty_         = p_x_->transpose() * (*p_y_);
      currentSet_.resize(p_x_->cols());
      for(int i = currentSet_.begin(); i < currentSet_.end(); i++)
        currentSet_[i] = i;
      // initialize penalty term
      p_penalty_->update(currentBeta_);
      return updateSolver();
  }
  /*Initialization of the solver */
  STK::Real LassoSolver::updateSolver()
  {
#ifdef HD_DEBUG
    std::cout << "Entering LassoSolver::updateSolver." << std::endl;
    //check the existence pointers to the data and the response
    if(p_x_ == 0)
      throw STK::invalid_argument(STK::String("p_x_ has not be set"));
    if(p_y_ == 0)
      throw STK::invalid_argument(STK::String("p_y_ has not be set"));
#endif
    b_.resize(currentSet_.range());
    for(int i = currentSet_.begin(); i < currentSet_.end(); i++)
    { b_[i]=Xty_[currentSet_[i]];}
    p_penalty_->update(currentBeta_);
    b_ *= p_penalty_->sqrtInvPenalty();
    x0_ = b_;
#ifdef HD_DEBUG
      std::cout << "LassoSolver::updateSolver done. llc =" << computeLlc() << std::endl;
#endif
      return computeLlc();

  }

  /*run the update of the penalty*/
  void LassoSolver::update(bool toUpdate)
  {
#ifdef HD_DEBUG
    if(p_penalty_ == 0)
      throw STK::invalid_argument(STK::String("p_penalty_ has not be set"));
#endif
    //if toUpdate check if some beta can be zeroed
    if (toUpdate) { updateCurrentBeta();}
    else
    {
      // otherwise just update current values of beta, penalty and b_
      for(int i = currentBeta_.begin(); i < currentBeta_.end(); i++)
      { p_beta_->elt(currentSet_[i]) = currentBeta_[i];}
      updateB();
    }
  }
  /* run the solver of the M step
   * @return
   */
  STK::Real LassoSolver::run(bool toUpdate)
  {
    //run the conjugate gradient
    int cgiter = cgsolver_.run();
#ifdef HD_DEBUG
      std::cout << "In LassoSolver::run. cgsolver_ run in " << cgiter << " iterations.\n";
#endif
    //backtransform the solution x to beta
    x0_ = cgsolver_.x();
    currentBeta_ = p_penalty_->sqrtInvPenalty() * x0_;
    //compute llc
    return computeLlc();
  }
  /* Computation of the completed loglikelihood*/
  STK::Real LassoSolver::computeLlc() const
  {
	  return(-(( *p_y_ - (currentX_ * currentBeta_)).norm2()/p_penalty_->sigma2() + p_penalty_->penaltyTerm(currentBeta_))/2.0);

  }

  /* compute an initial value of beta using ols */
  void LassoSolver::computeInitialBeta()
  {
    // compute initial value using OLS
    Xty_ = currentX_.transpose() * *p_y_;
    InitLassoMultiplicator mult(p_x_, p_penalty_->lambda());
    STK::CG<InitLassoMultiplicator,STK::VectorX,InitFunctor> cginitial(mult, Xty_, 0, 1e-05);
    cginitial.run();
    p_beta_->move(cginitial.x());
    disruptsBeta();
    currentBeta_ = *p_beta_;
  }

  /* update all the current variables*/
  void LassoSolver::updateCurrentBeta()
  {
#ifdef HD_DEBUG
      std::cout << "Entering LassoSolver::updateCurrentBeta." << std::endl;
      std::cout << "currentSet_.range() =" << currentSet_.range() << std::endl;
#endif
    int nbActiveVariables = currentBeta_.size();
    //thresholding of the new estimates, currentbeta become 0.
    for(int i = currentBeta_.begin(); i < currentBeta_.end(); i++)
    {
      if (std::abs(currentBeta_[i]) < threshold_)
      { currentBeta_[i] = 0.; nbActiveVariables--;
#ifdef HD_VERBOSE_DEBUG
        std::cout << "Removing variable " << currentSet_[i] << std::endl;
#endif
      }
      // copy current value to beta
      p_beta_->elt(currentSet_[i]) = currentBeta_[i];
    }
    //if TRUE, at least one variable became 0 in the thresholding step
    if(nbActiveVariables != currentBeta_.size())
    {
#ifdef HD_DEBUG
      std::cout << "In updateCurrentBeta. Removing " << currentBeta_.size() - nbActiveVariables << " variables\n";
#endif
      //erase index of the zero estimates
      for(int i = currentBeta_.lastIdx(); i >= currentBeta_.begin(); i--)
      { if (currentBeta_[i]==0.) { currentSet_.erase(i); x0_.erase(i);}}
#ifdef HD_DEBUG
      std::cout << "In updateCurrentBeta. currentSet_.range() =" << currentSet_.range() << std::endl;
#endif
      // update currentX, currentBeta, penalty and b_
      updateSystem();
    }
    else // just update penalty and b_
    { updateB();}
  }

  /* update the second member */
  void LassoSolver::updateB()
  {
    //compute the b of the linear system Ax=b
    //b_ = p_penalty_->sqrtInvPenalty() * currentX_.transpose() * (*p_y_);
    b_.resize(currentSet_.range());
    for(int i = currentSet_.begin(); i < currentSet_.end(); i++)
    { b_[i]=Xty_[currentSet_[i]];}
    p_penalty_->update(currentBeta_);
    b_ *= p_penalty_->sqrtInvPenalty();
  }

  /* update the second member */
  void LassoSolver::updateSystem()
  {
    b_.resize(currentSet_.range());
    currentBeta_.resize(currentSet_.range());
    currentX_.resize(p_x_->rows(),currentSet_.range());
    for(int i = currentSet_.begin(); i < currentSet_.end(); i++)
    {
      b_[i]            = Xty_[currentSet_[i]];
      currentBeta_[i]  = p_beta_->elt(currentSet_[i]);
      currentX_.col(i) = p_x_->col(currentSet_[i]);
    }
    p_penalty_->update(currentBeta_);
    b_ *= p_penalty_->sqrtInvPenalty();
#ifdef HD_VERBOSE_DEBUG
      std::cout << "In updateCurrentBeta. x0_.range() =" << x0_.range() << std::endl;
      std::cout << "In updateSystem. b_.range() =" << b_.range() << std::endl;
      std::cout << "In updateSystem. currentBeta_.range() =" << currentBeta_.range() << std::endl;
      std::cout << "In updateSystem. currentX_.rows() =" << currentX_.rows() << std::endl;
      std::cout << "In updateSystem. currentX_.cols() =" << currentX_.cols() << std::endl;
#endif
  }
} // HD
