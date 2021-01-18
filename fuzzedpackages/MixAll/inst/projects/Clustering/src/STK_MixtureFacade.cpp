/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 4 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_MixtureFacade.cpp
 *  @brief In this file we implement the StrategyFacade class.
 **/

#include "../include/STK_IMixtureComposer.h"
#include "../include/STK_IMixtureLearner.h"
#include "../include/STK_MixtureFacade.h"
#include "../include/MixtureInit/STK_MixtureInit.h"
#include "../include/MixtureStrategy/STK_MixtureStrategy.h"
#include "../include/MixtureAlgo/STK_MixtureAlgo.h"
#include "../include/MixtureAlgo/STK_MixtureAlgoLearn.h"

namespace STK
{

StrategyFacade::~StrategyFacade() { if (p_strategy_) delete p_strategy_;}

void StrategyFacade::createSimpleStrategy( Clust::initType init, int nbTryInInit, Clust::algoType initAlgo, int nbInitIter, Real initEpsilon
                                        , int nbTry, Clust::algoType algo, int nbIter, Real epsilon)
{
  IMixtureInit* p_init = Clust::createInit(init, nbTryInInit, initAlgo, nbInitIter, initEpsilon);
  IMixtureAlgo* p_algo = Clust::createAlgo(algo, nbIter, epsilon);

  p_strategy_ = Clust::createSimpleStrategy(p_model_, nbTry, p_init, p_algo);
}

/* create a FullStrategy */
void StrategyFacade::createFullStrategy( Clust::initType init, int nbTryInInit, Clust::algoType initAlgo, int nbInitIter, Real initEpsilon
                       , int nbTry, int nbInitRun, int nbShortRun
                       , Clust::algoType shortAlgo, int nbShortIter, Real shortEpsilon
                       , Clust::algoType longAlgo, int nblongIter, Real longEpsilon)
{
  IMixtureInit* p_init = Clust::createInit(init, nbTryInInit, initAlgo, nbInitIter, initEpsilon);
  IMixtureAlgo* p_shortAlgo = Clust::createAlgo(shortAlgo, nbShortIter, shortEpsilon);
  IMixtureAlgo* p_longAlgo = Clust::createAlgo(longAlgo, nblongIter, longEpsilon);

  p_strategy_ = Clust::createFullStrategy(p_model_, nbTry, nbInitRun, p_init, nbShortRun, p_shortAlgo, p_longAlgo);
}

bool StrategyFacade::run()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("------------------------------\n")
           << _T("Entering StrategyFacade::run()\n");
#endif
  bool flag = false;
  if (p_strategy_)
  {
    if (p_strategy_->run()) { flag = true;}
    else
    {
      msg_error_ += STKERROR_NO_ARG(StrategyFacade::run,strategy failed\n);
      msg_error_ += p_strategy_->error();
#ifdef STK_MIXTURE_VERBOSE
      stk_cout << _T("StrategyFacade:run() terminated without success.\n")
               << msg_error_ << _T("\n")
               << _T("------------------------------------------------\n");
#endif
    }
    // p_model_->imputationStep();
    p_model_->finalizeStep();
  }
  else
  { msg_error_ = STKERROR_NO_ARG(MixtureFacade::run(),strategy is not created);}
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("StrategyFacade:run() terminated.\np_model->lnLikelihood() =")
           << p_model_->lnLikelihood() << _T("\n")
           << _T("--------------------------------\n");
#endif
  return flag;
}

LearnFacade::~LearnFacade() { if (p_algo_) delete p_algo_;}

void LearnFacade::createImputeAlgo(  Clust::algoLearnType algo, int nbIter, Real epsilon)
{
  p_algo_ = Clust::createLearnAlgo(algo, nbIter, epsilon);
  p_algo_->setModel(p_model_);
}

/* create a FullStrategy */
void LearnFacade::createSimulAlgo( Clust::algoLearnType algo, int nbIter)
{
  p_algo_ = Clust::createLearnAlgo(algo, nbIter, 0.);
  p_algo_->setModel(p_model_);
}

bool LearnFacade::run()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("------------------------------\n")
           << _T("Entering LearnFacade::run()   \n");
#endif
  // initialize model
  p_model_->initializeStep();
  // start algorithm
  bool flag = false;
  if (p_algo_)
  {
    if (p_algo_->run()) { flag = true;}
    else
    {
      msg_error_ += STKERROR_NO_ARG(StrategyFacade::run,strategy failed\n);
      msg_error_ += p_algo_->error();
#ifdef STK_MIXTURE_VERBOSE
      stk_cout << _T("LearnFacade:run() terminated without success.   \n")
               << msg_error_ << _T("\n")
               << _T("------------------------------------------------\n");
#endif
    }
    // finalize any
    p_model_->finalizeStep();
  }
  else
  {
    msg_error_ = STKERROR_NO_ARG(MixtureFacade::run(),algo is not created);
  }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("StrategyFacade:run() terminated.\n")
           << _T("p_model->lnLikelihood() =") << p_model_->lnLikelihood() << _T("\n")
           << _T("--------------------------------\n");
#endif
  return flag;
}

}  // namespace STK




