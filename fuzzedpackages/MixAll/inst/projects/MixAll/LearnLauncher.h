/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2015  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 10 May 2016
 * Author:   Iovleff, serge.iovleff@stkpp.org
 **/

/** @file LearnLauncher.h
 *  @brief In this file we define the LearnLauncher class which
 *  construct properly a mixture model.
 **/


#ifndef STK_LEARNLAUNCHER_H
#define STK_LEARNLAUNCHER_H

#include "ILauncher.h"

namespace STK
{

/**The LearnLauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
class LearnLauncher: public ILauncher
{
  public:
    /** constructor.
     * @param model a reference on the current model
     * @param algo the algorithm defined in R
     **/
    LearnLauncher( Rcpp::S4 model, Rcpp::CharacterVector models, Rcpp::S4 algo );
    /** constructor with a list of component.
     *  @param model a reference on the current model
     *  @param algo the algorithm defined in R
     **/
    LearnLauncher( Rcpp::S4 model, Rcpp::S4 algo );
    /** destructor. */
    virtual ~LearnLauncher();
    /** run the estimation */
    bool run();

  protected:
    /** strategy from the R side */
    Rcpp::S4    s4_algo_;
    /** character string with the model selection criterion name */
    String criterion_;
    /** learning algorithm to run */
    IMixtureAlgoLearn* p_algo_;
    /** criterion to run */
    IMixtureCriterion* p_criterion_;

  private:
    /** pointer on the main learner */
    IMixtureLearner* p_learner_;
    /** Is the model with mixed data ? */
    bool isMixedData_;

    /** Select the best model among the models and nbCluster given.
     *  @return the value of the best criteria.
     **/
    Real selectBestSingleModel();
    /** Select the best model among the models and nbCluster given.
     *  @return the value of the best criteria.
     **/
    Real selectBestMixedModel();
};

} // namespace STK

#endif /* STK_LEARNLAUNCHER_H */
