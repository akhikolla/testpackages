/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff, University Lille 1, Inria

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 4 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_MixtureFacade.cpp
 *  @brief In this file we implement the ClusterFacade class.
 **/

#include <MixAll/ClusterFacade.h>

namespace STK
{

ClusterFacade::~ClusterFacade() { if (p_strategy_) delete p_strategy_;}


/* create a FullStrategy */
void ClusterFacade::createFullStrategy(Rcpp::S4 s4_strategy)
{
  // get fields of the strategies
  int nbTry = s4_strategy.slot("nbTry");
  int nbShortRun = s4_strategy.slot("nbShortRun");
  Rcpp::S4 R_initMethod = s4_strategy.slot("initMethod");
  Rcpp::S4 R_shortAlgo = s4_strategy.slot("shortAlgo");
  Rcpp::S4 R_longAlgo = s4_strategy.slot("longAlgo");

  // get fields of the initMethod
  String method = R_initMethod.slot("method");
  STK::Clust::initType init = STK::Clust::stringToInit(method);
  int nbInit = R_initMethod.slot("nbInit");
  Rcpp::S4 R_initAlgo = R_initMethod.slot("algo");

  // get fields of the initAlgo
  String initAlgoName = R_initAlgo.slot("algo");
  STK::Clust::algoType initAlgo = STK::Clust::stringToAlgo(initAlgoName);
  int nbInitIter = R_initAlgo.slot("nbIteration");
  STK::Real initEpsilon = R_initAlgo.slot("epsilon");

  // get fields of the shortAlgo
  String shortAlgoName = R_shortAlgo.slot("algo");
  STK::Clust::algoType shortAlgo = STK::Clust::stringToAlgo(shortAlgoName);
  int nbShortIter = R_shortAlgo.slot("nbIteration");
  STK::Real shortEpsilon = R_shortAlgo.slot("epsilon");

  // get fields of the longAlgo
  String longAlgoName = R_longAlgo.slot("algo");
  STK::Clust::algoType longAlgo = STK::Clust::stringToAlgo(longAlgoName);
  int nbLongIter = R_longAlgo.slot("nbIteration");
  STK::Real longEpsilon = R_longAlgo.slot("epsilon");

  // create STK objects
  STK::IMixtureInit* p_init = STK::Clust::createInit(init, 1, initAlgo, nbInitIter, initEpsilon);
  STK::IMixtureAlgo* p_shortAlgo = STK::Clust::createAlgo(shortAlgo, nbShortIter, shortEpsilon);
  STK::IMixtureAlgo* p_longAlgo = STK::Clust::createAlgo(longAlgo, nbLongIter, longEpsilon);
  p_strategy_ = STK::Clust::createFullStrategy(p_model_, nbTry, nbInit, p_init, nbShortRun, p_shortAlgo, p_longAlgo);
}

bool ClusterFacade::run()
{
  bool flag = false;
  if (p_strategy_)
  {
    flag = p_strategy_->run();
    if (!flag) { msg_error_ = p_strategy_->error();}
    p_model_->imputationStep(); // impute the missing values
    p_model_->finalizeStep();   // finalize any model stuff needed
  }
  else
  { msg_error_ = STKERROR_NO_ARG(ClusterFacade::run,strategy is not set);}
  return flag;
}

} // namespace STK

