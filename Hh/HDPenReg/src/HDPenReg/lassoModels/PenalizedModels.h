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

/** @file PenalizedModels.h
 *  @brief In this file, we find the definition of the PenalizedModels class.
 **/


#ifndef PENALIZEDMODELS_H_
#define PENALIZEDMODELS_H_

#include <RTKpp.h>
#include <StatModels.h>
#include "IPenalizedSolver.h"
#include "IPenalty.h"

namespace HD
{

  template <class Model> struct ModelTraits;

  template<class Model>
  class PenalizedModels : public STK::IStatModelBase
  {
    public:

    typedef typename ModelTraits<Model>::Solver Solver;
    typedef typename ModelTraits<Model>::Penalty Penalty;

    /**
     * Constructor
     * @param p_x each row contains values of covariates for an individuals
     * @param y response of each individuals
     * @param beta Initialization of beta_ (least-Square solution)
     * @param p_penalty pointer on the penalty of the model
     * @param p_solver pointer on the solver to use in the Mstep
     * @param threshold tolerance for thresholding the solution
     **/
    PenalizedModels( STK::ArrayXX const* p_x, STK::VectorX const* p_y
                   , STK::VectorX const& beta
                   , Solver* p_solver = 0)
                   : STK::IStatModelBase(p_x->sizeRows(),p_x->sizeCols())
                   , p_x_(p_x)
                   , p_y_(p_y)
                   , beta_(beta)
                   , p_solver_(p_solver)
                   , p_penalty_(0)
    {}

    PenalizedModels( STK::ArrayXX const* p_x, STK::VectorX const* p_y)
                   : STK::IStatModelBase(p_x->sizeRows(),p_x->sizeCols())
                   , p_x_(p_x)
                   , p_y_(p_y)
                   , beta_(p_x->sizeCols())
                   , p_solver_(0)
                   , p_penalty_(0)
    {}

    PenalizedModels(): STK::IStatModelBase(),p_x_(0),p_y_(0),beta_(),p_solver_(0), p_penalty_(0)
    {}

    /** destructor*/
    ~PenalizedModels()
    {
       if(p_penalty_) delete p_penalty_;
       if(p_solver_) delete p_solver_;
    };
    //getter
    /**@return a pointer to the data*/
    inline STK::ArrayXX const* p_x() const {return p_x_;}
    /** @return a pointer to the response */
    inline STK::VectorX const* p_y() const {return p_y_;}
    /**@return the estimated beta*/
    inline STK::VectorX const& beta() const {return beta_;}
    /**@return the estimated beta*/
    inline STK::Real beta(int i) const {return beta_[i];}
    /**@return p_solver_*/
    inline Solver const* p_solver() const {return p_solver_;}
    /**@return p_penalty_*/
    inline Penalty const* p_penalty() const {return p_penalty_;}

    // short cuts for the temporaries arrays used in solver
    /** @return  the non zero beta*/
    inline STK::VectorX const& currentBeta() const {return p_solver_->currentBeta();}
    /** @return  the non zero beta*/
    inline STK::VectorX const* p_currentBeta() const {return p_solver_->p_currentBeta();}
    /** @return  data*/
    inline STK::ArrayXX const& currentX() const {return p_solver_->currentX();}
    /** @return  data*/
    inline STK::ArrayXX const* p_currentX() const {return p_solver_->p_currentX();}
    /**@return currentSet_*/
    inline STK::VectorXi const& currentSet() const {return p_solver_->currentSet();};

    //setter
    inline void setP_y(STK::VectorX const* p_y) {p_y_ = p_y;}
    inline void setP_x(STK::ArrayXX const* p_x)
    {
      p_x_ = p_x;
      this->setNbSample(p_x_->sizeRows());
      this->setNbVariable(p_x_->sizeCols());
    }
    //setter
    /** set the solver
     * @param p_solver a pointer to a IPenalizedSolver object
     */
    inline void setSolver(Solver* p_solver) { p_solver_=p_solver;}
    /** EStep update the current beta if toUpdate is @c true
     *  and update the latent variable in the penalty term*/
    inline void eStep(bool toUpdate)
    {
#ifdef HD_DEBUG
      if(p_solver_ == 0)
      { throw STK::invalid_argument(STK::String("p_solver_ has not be set"));}
#endif
      p_solver_->update(toUpdate);
    }
    /** MStep update the estimation of beta*/
    inline void mStep(bool toUpdate)
    {
#ifdef HD_DEBUG
      if(p_solver_ == 0)
      { throw STK::invalid_argument(STK::String("p_solver_ has not be set"));}
#endif
      this->setLnLikelihood(p_solver_->run(toUpdate));
    }
    /** update solver and set the initial beta0 */
    STK::Real updateSolver()  { return p_solver_->updateSolver();}
    /**initialize the solver */
    STK::Real initializeSolver()  { return p_solver_->initializeSolver();}

    protected:
      ///pointer to the data
      STK::ArrayXX const* p_x_;
      ///response (length n_)
      STK::VectorX const* p_y_;
      ///estimates of the model (length p_)
      STK::VectorX beta_;
      ///pointer to the solver for solving the M step
      Solver* p_solver_;
      ///pointer to the penalty of the model
      Penalty* p_penalty_;
  };
}

#endif /* PENALIZEDMODELS_H_ */
