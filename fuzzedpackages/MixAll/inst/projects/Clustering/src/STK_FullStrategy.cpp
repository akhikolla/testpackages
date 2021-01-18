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
 * created on: 3 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_FullStrategy.cpp
 *  @brief In this file we implement the FullStrategy class
 **/

#include <Sdk.h>

#include "../include/MixtureStrategy/STK_FullStrategy.h"
#include "../include/MixtureInit/STK_MixtureInit.h"
#include "../include/MixtureAlgo/STK_MixtureAlgo.h"
#include "../include/STK_IMixtureComposer.h"

namespace STK
{

/* destructor */
FullStrategyParam::~FullStrategyParam()
{
  if (p_shortAlgo_) delete p_shortAlgo_;
  if (p_longAlgo_) delete p_longAlgo_;
}


/* run the full strategy */
bool FullStrategy::run()
{
  // add some perturbation to the tik and compute the ln-likelihood
  p_model_->setState(Clust::modelInitialized_);
#ifdef STK_MIXTURE_DEBUG
  p_model_->writeParameters(stk_cout);
#endif
  p_model_->randomFuzzyInit();
#ifdef STK_MIXTURE_DEBUG
  p_model_->writeParameters(stk_cout);
#endif
  Real initialValue = p_model_->lnLikelihood();
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("<+++\n");
  stk_cout << _T("Entering FullStrategy::run() with nbTry_ = ") << nbTry_
           << _T(", nbShortRun_ = ") << p_param_->nbShortRun_
           << _T(", p_model_->lnLikelihood() = ") << initialValue
           << _T("\n");
#endif
  IMixtureComposer *p_bestModel =0, *p_bestShortModel =0;
  // start estimation
  try
  {
    // Main loop. If the Full strategy success in estimating a model, the
    // iterations are stopped and the best model find is stored in p_model_
    for (int iTry = 0; iTry < nbTry_; ++iTry)
    {
      // in case nbShortRun_==0: initialize directly p_bestShortModel
      if (p_param_->nbShortRun_ <= 0)
      {
        if (!initStep(p_bestShortModel))
        {
          msg_error_ += STKERROR_NO_ARG(FullStrategy::run,init step failed\n);
          msg_error_ += p_param_->p_shortAlgo_->error();
#ifdef STK_MIXTURE_VERBOSE
          stk_cout << _T("In FullStrategy::run()") << _T(", iTyry =") << iTry
                   << _T(", init step failed\n");
          stk_cout << msg_error_ << _T("\n");
#endif
        }
      }
      else
      {
#ifdef STK_MIXTURE_VERY_VERBOSE
        stk_cout << _T("In FullStrategy::run(), entering short run steps\n")
                 << _T("iTyry =") << iTry << _T("\n");
#endif
        Real valueBest = -Arithmetic<Real>::infinity();
        for (int iShort=0; iShort < p_param_->nbShortRun_; ++iShort)
        {
#ifdef STK_MIXTURE_VERY_VERBOSE
          stk_cout << _T("In FullStrategy::run(), iShort =") << iShort << _T("\n");
#endif
          // perform nbInitRun_ initialization step and get the best result in p_bestModel
          if (!initStep(p_bestModel))
          {
            msg_error_ += STKERROR_NO_ARG(FullStrategy::run,init step failed\n);
            msg_error_ += p_param_->p_shortAlgo_->error();
#ifdef STK_MIXTURE_VERBOSE
            stk_cout << _T("In FullStrategy::run()") << _T(", iTyry =") << iTry << _T(", iShort =") << iShort
                     << _T(", init step failed\n");
            stk_cout << msg_error_ << _T("\n");
#endif
          }
          // if we get a better result, store it in p_bestShortModel
          Real value = p_bestModel->lnLikelihood();
          if( valueBest<value)
          {
            std::swap(p_bestShortModel, p_bestModel);
            valueBest  = value;
#ifdef STK_MIXTURE_VERY_VERBOSE
            stk_cout << _T("In FullStrategy::run()")
                     << _T(", iTyry =") << iTry << _T(", iShort =") << iShort
                     << _T(", get better value in short run. valueBest =") << valueBest << _T("\n");
#endif
          }
        } // ishort
        // release memory
        if (p_bestModel) { delete p_bestModel; p_bestModel = 0;}
      }
      // in case all initialization failed
      if (!p_bestShortModel) { p_bestShortModel = p_model_->clone();}
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("In FullStrategy::run() all short run done") << _T(", iTyry =") << iTry  << _T(" terminated.\n")
           << _T("p_bestShortModel->lnLikelihood() = ") << p_bestShortModel->lnLikelihood()
           << _T("\n");
#endif
      // start a long run with p_bestShortModel. If success, save model
      // and exit the iTry loop
      p_param_->p_longAlgo_->setModel(p_bestShortModel);
      if (!p_param_->p_longAlgo_->run())
      {
        msg_error_ += STKERROR_NO_ARG(FullStrategy::run,long algo failed\n);
        msg_error_ += p_param_->p_longAlgo_->error();
#ifdef STK_MIXTURE_VERBOSE
        stk_cout << _T("In FullStrategy::run(): Long Algo failed\n");
#endif
      }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("In FullStrategy::run() long run") << _T(", iTyry =") << iTry  << _T(" terminated.\n")
           << _T("p_bestShortModel->lnLikelihood() = ") << p_bestShortModel->lnLikelihood()
           << _T("\n");
#endif
      // if we get a better result, store it in p_model_ and stop to try
      if( p_model_->lnLikelihood()<p_bestShortModel->lnLikelihood())
      {
        std::swap(p_model_, p_bestShortModel);
        break;
      }
#ifdef STK_MIXTURE_VERBOSE
      stk_cout << _T("In FullStrategy::run(), iTry =") << iTry << _T(" failed\n");
#endif
      // release memory before next try
      if (p_bestModel)      delete p_bestModel;
      p_bestModel = 0;
      if (p_bestShortModel) delete p_bestShortModel;
      p_bestShortModel = 0;
    } // end iTry
  }
  catch (Exception const& e)
  {
    if (p_bestModel)      delete p_bestModel;
    if (p_bestShortModel) delete p_bestShortModel;
    msg_error_ += e.error();
    return false;
  }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << "FullStrategy::run() terminated. \n";
  stk_cout << _T("+++>\n");
#endif
  // normally not needed
  if (p_bestModel) delete p_bestModel;
  if (p_bestShortModel)   delete p_bestShortModel;
  if (p_model_->lnLikelihood() <= initialValue)
  {
    msg_error_ += STKERROR_NO_ARG(FullStrategy::run,No gain\n);
    return false;
  }
  return true;
}

/* Perform the Initialization step*/
bool FullStrategy::initStep(IMixtureComposer*& p_bestModel)
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("<+++++\n");
  stk_cout << _T("Entering FullStrategy::initStep\n");
  stk_cout << _T("nbInitRun = ") <<  p_param_->nbInitRun_ << _T("\n");
#endif
  bool flag = true;
  IMixtureComposer* p_initModel = 0;
  try
  {
    Real valueBest = -Arithmetic<Real>::infinity();
    for (int iInitRun=0; iInitRun < p_param_->nbInitRun_; iInitRun++)
    {
      //  Initialize a new model if necessary
      if (!p_initModel)
      { p_initModel = p_model_->create();}
      p_init_->setModel(p_initModel);
      if (!p_init_->run())
      {
#ifdef STK_MIXTURE_VERBOSE
        stk_cout << _T("FullStrategy::initStep, iInitRun=") << iInitRun
                 << _T(", initialization failed:\n");
        stk_cout << p_init_->error() << _T("\n");
#endif
        msg_error_ += p_init_->error();
      }
      else
      {
        // if we get a better result, swap initModel with bestModel
        Real value = p_initModel->lnLikelihood();
        if( (valueBest < value) && isFinite(value))
        {
          std::swap(p_initModel, p_bestModel);
          valueBest = value;
#ifdef STK_MIXTURE_VERY_VERBOSE
        stk_cout << _T("FullStrategy::initStep, iInitRun =") << iInitRun
                 << _T(", get a better model with value =") << valueBest << _T("\n");
#endif
        }
      }
    } // iInitRun
    // In case we never get a better model, clone current model
    // and perform short run with the current model
    if (!p_bestModel) { p_bestModel = p_model_->clone();}
    p_param_->p_shortAlgo_->setModel(p_bestModel);
    if (!p_param_->p_shortAlgo_->run())
    {
      msg_error_ = STKERROR_NO_ARG(FullStrategy::initStep,short algo failed\n);
      msg_error_ = p_param_->p_shortAlgo_->error();
#ifdef STK_MIXTURE_VERBOSE
      stk_cout << _T("In FullStrategy::initStep() shortAlgo failed:\n");
      stk_cout << msg_error_ << _T("\n");
#endif
    }
  }
  catch (Exception const& e)
  {
    msg_error_ = e.error();
    flag = false;
  }
  // in case all initialization failed or nbInitRun_ <= 0
  if (!p_bestModel) { p_bestModel = p_model_->clone();}
  if (p_initModel)
  {
#ifdef STK_MIXTURE_DEBUG_CREATE
        stk_cout << _T("FullStrategy::initStep terminated. Deleting p_initModel.\n");
#endif
    delete p_initModel; p_initModel = 0;
#ifdef STK_MIXTURE_DEBUG_CREATE
        stk_cout << _T("FullStrategy::initStep p_initModel deleted.\n");
        stk_cout << _T("p_bestModel->writeParameters\n");
        p_bestModel->writeParameters(stk_cout);
#endif
  }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("FullStrategy::initStep done\n");
  stk_cout << _T("p_bestModel->lnLikelihood() = ") <<  p_bestModel->lnLikelihood() << _T("\n");
  stk_cout << _T("+++++>\n");
#endif
  return flag;
}

} // namespace STK



