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

/** @file LogisticLassoSolver.cpp
 *  @brief In this file .
 **/


#include "LogisticLassoSolver.h"

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
  LogisticLassoSolver::LogisticLassoSolver( STK::ArrayXX const* p_x, STK::VectorX const* p_y, STK::VectorX* p_beta
                                          , STK::Real const& threshold, STK::Real const& epsCG
                                          , LassoPenalty* p_penalty)
                                          : IPenalizedSolver(p_beta, p_x, p_y, threshold)
                                          , z_(), b_(), x0_()
                                          , p_penalty_(p_penalty)
                                          , mult_(), cgsolver_(), cginit_()
  {
    // currentX_  is initialized to *p_x in IPenalizedSolver
    // currentSet_ is initialized in IPenalizedSolver to currentSet_[i] = i
    z_ = *p_y_;
    // currentBeta_ is initialized in computeInitialBeta
    computeInitialBeta();
#ifdef HD_DEBUG
    std::cout << "In LogisticLassoSolver::LogisticLassoSolver initial beta computed\n";
#endif
    // initialize CG
    cgsolver_.setMultFunctor(mult_);
    cgsolver_.setEps(epsCG);
    cgsolver_.setB(b_);
    cginit_.p_x0_ = &x0_;
    cgsolver_.setInitFunctor(&cginit_);
    //intialize mult functor for CG
    mult_.p_x_          = p_currentX();
    mult_.p_sigma2_     = p_penalty_->p_sigma2();
    mult_.p_sqrtInvPenalty_ = p_penalty_->p_sqrtInvPenalty();
    // initialize penalty term
    p_penalty_->update(currentBeta_);
    x0_ = p_penalty_->sqrtInvPenalty() * (currentX_.transpose() * z_);
#ifdef HD_DEBUG
    std::cout << "LogisticLassoSolver Initialized. Likelihood =" << computeLlc();
#endif
  }

  /*Initialization of the solver*/
  STK::Real LogisticLassoSolver::initializeSolver()
  {
#ifdef HD_DEBUG
    //check the existence pointers to the data and the response
    if(p_x_ == 0)
      throw STK::invalid_argument(STK::String("p_x_ has not be set"));
    if(p_y_ == 0)
      throw STK::invalid_argument(STK::String("p_y_ has not be set"));
#endif
    z_ = *p_y_;
    return updateSolver();
  }
  /* initialize beta0 */
  STK::Real LogisticLassoSolver::updateSolver()
  {
    // reinitialize all containers
    currentX_    = *p_x_;
    currentBeta_ = *p_beta_;
    currentSet_.resize(p_x_->cols());
    for(int i = currentSet_.begin(); i < currentSet_.end(); i++) currentSet_[i] = i;
    // initialize penalty term
    p_penalty_->update(currentBeta_);
    updateZ();
    x0_  = p_penalty_->sqrtInvPenalty() * (currentX_.transpose() * z_);
    return computeLlc();
  }


  /* Computation of the completed loglikelihood*/
  STK::Real LogisticLassoSolver::computeLlc() const
  {
    return -( ((z_ - (currentX_ * currentBeta_) ).square().sum())/p_penalty_->sigma2()
              + p_penalty_->penaltyTerm(currentBeta_)
            )/2;;
  }

  /*run the update of the penalty*/
  void LogisticLassoSolver::update(bool toUpdate)
  {
#ifdef HD_DEBUG
    if(p_penalty_ == 0)
      throw STK::invalid_argument(STK::String("p_penalty_ has not be set"));
#endif
    //check if we can reduce the dimension
    if (toUpdate)
    { updateCurrentBeta();}
    else
    {
      // otherwise just update current values of beta
      for(int i = currentBeta_.begin(); i < currentBeta_.end(); i++)
      { p_beta_->elt(currentSet_[i]) = currentBeta_[i];}
    }
    // update penalty
    p_penalty_->update(currentBeta_);
    // update z_ values
    updateZ();
  }
  /* run the solver of the M step
   * @return
   */
  STK::Real LogisticLassoSolver::run(bool toUpdate)
  {
    //update the b of the linear system Ax=b and run cg
    b_ = p_penalty_->sqrtInvPenalty() * (currentX_.transpose() * z_);
    int cgiter = cgsolver_.run();
#ifdef HD_DEBUG
      std::cout << "In LogisticLassoSolver::run. cgsolver_ run in " << cgiter << " iterations.\n";
#endif
    //back-transform the solution x to beta
    x0_ = cgsolver_.x();
    currentBeta_ = p_penalty_->sqrtInvPenalty() * x0_;
    //compute llc
    return computeLlc();
  }

  /* compute an initial value of beta using ols */
  void LogisticLassoSolver::computeInitialBeta()
  {
    // compute initial value using OLS
    STK::Vector Xty = currentX_.transpose() * z_;
    InitLassoMultiplicator mult(p_x_, p_penalty_->lambda());
    STK::CG<InitLassoMultiplicator,STK::VectorX,InitFunctor> cginit(mult, Xty, 0, 1e-05);
    cginit.run();
    p_beta_->move(cginit.x());
    disruptsBeta();
    currentBeta_ = *p_beta_;
  }

  /* update all the current variables*/
  void LogisticLassoSolver::updateCurrentBeta()
  {
#ifdef HD_DEBUG
      std::cout << "Entering LogisticLassoSolver::updateCurrentBeta." << std::endl;
      std::cout << "threshold =" << threshold_ << std::endl;
      std::cout << "currentSet_.range() =" << currentSet_.range() << std::endl;
#endif
    int nbActiveVariables = currentBeta_.size();
    //thresholding of the new estimates
    for(int i = currentBeta_.begin(); i < currentBeta_.end(); i++)
    {
      if(std::abs(currentBeta_[i]) < threshold_)
      { currentBeta_[i] = 0; nbActiveVariables--;}
      // copy current value to beta
      p_beta_->elt(currentSet_[i]) = currentBeta_[i];
    }
    //if TRUE, at least one variable became 0 in the previous M step
    if(nbActiveVariables != currentBeta_.size())
    {
#ifdef HD_DEBUG
      std::cout << "In updateCurrentBeta. Removing " << currentBeta_.size() - nbActiveVariables << " variables\n";
//      std::cout << "x0_.range() =" << x0_.range() << std::endl;
#endif
      //erase index of the zero estimates
      for(int i = currentBeta_.lastIdx(); i >= currentBeta_.begin(); i--)
      { if (currentBeta_[i]==0.) { currentSet_.erase(i); x0_.erase(i);}}
      // update currentX and currentBeta
      updateSystem();
#ifdef HD_DEBUG
      std::cout << "In updateCurrentBeta removing terminated. currentSet_.range() =" << currentSet_.range() << std::endl;
#endif
    }
    // update penalty term
    p_penalty_->update(currentBeta_);
  }

  /* update the system */
  void LogisticLassoSolver::updateSystem()
  {
    //resize currentBeta_ and currentX_ with only the active variable
    currentBeta_.resize(currentSet_.range());
    currentX_.resize(p_x_->rows(), currentSet_.range());
    for(int i = currentBeta_.begin(); i < currentBeta_.end(); i++)
    {
      currentBeta_[i]  = p_beta_->elt(currentSet_[i]);
      currentX_.col(i) = p_x_->col(currentSet_[i]);
    }
  }

  /* update z values */
  void LogisticLassoSolver::updateZ()
  {
    for(int i = z_.begin(); i < z_.end(); i++)
    {
      STK::Real aux = currentBeta_.dot(currentX_.row(i));
      STK::Real det = ((*p_y_)[i] == 1) ? (normal_.cdf(aux)) : -normal_.cdf(-aux);
      if (std::abs(det) < 1e-10) det = 1e-10;
      z_[i] = aux +  normal_.pdf(aux)/ det;
    }
  }


} // namespace HD

