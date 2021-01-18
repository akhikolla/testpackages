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

/** @file STK_XemStrategy.cpp
 *  @brief In this file we implement the XemStrategy class
 **/

#include <Sdk.h>

#include "../include/MixtureStrategy/STK_XemStrategy.h"
#include "../include/MixtureInit/STK_MixtureInit.h"
#include "../include/MixtureAlgo/STK_MixtureAlgo.h"
#include "../include/STK_IMixtureComposer.h"

namespace STK
{


/* destructor */
XemStrategyParam::~XemStrategyParam()
{
  if (p_shortAlgo_) delete p_shortAlgo_;
  if (p_longAlgo_) delete p_longAlgo_;
}


/* run the xem strategy */
bool XemStrategy::run()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Entering XemStrategy::run() with:\n")
           << _T("nbTry_ = ") << nbTry_ << _T("\n")
           << _T("nbShortRun_ = ") << p_param_->nbShortRun_ << _T("\n");
#endif
  // the current model is used in the short runs
  IMixtureComposer* p_currentModel     = p_model_->create();
  IMixtureComposer* p_currentBestModel = p_model_->create();
  // add some perturbation to the tik and compute the ln-likelihood
  if (p_model_->state() < 1) { p_model_->randomFuzzyInit();}
  Real value = p_model_->lnLikelihood();
  // start estimation
  try
  {
    for (int iTry = 0; iTry < nbTry_; ++iTry)
    {
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("-------------------------------\n")
           << _T("try number = ") << iTry << _T("\n");
#endif
      // find best of the shortModel and save it in p_currentBestModel
      for (int iShortRun = 0; iShortRun < p_param_->nbShortRun_; ++iShortRun)
      {
        // initialize current model
        p_init_->setModel(p_currentModel);
        if (p_init_->run())
        {
          // perform short run on the current model
          p_param_->p_shortAlgo_->setModel(p_currentModel);
          p_param_->p_shortAlgo_->run();
          // if we get a better result, swap it with currentBestModel
          if( p_currentBestModel->lnLikelihood()<p_currentModel->lnLikelihood())
          { std::swap(p_currentModel, p_currentBestModel);}
        } // initialization
      } // iShortRun
      // in case nbShortRun_==0
      // try to initialize bestCurrentModel, otherwise go to a next try
      if (p_param_->nbShortRun_ == 0)
      {
        // initialize current model
        p_init_->setModel(p_currentBestModel);
        if (!p_init_->run())
        { continue; }// model not initialized, we go to the next trial
      }
#ifdef STK_MIXTURE_VERBOSE
      stk_cout << _T("iTry =") << iTry
               << _T(". In XemStrategy::run(), short run terminated. best model:\n");
      p_currentBestModel->writeParameters(stk_cout);
      stk_cout << _T("\n\n");
#endif
      // start a long run with the better model
      p_param_->p_longAlgo_->setModel(p_currentBestModel);
      if (p_param_->p_longAlgo_->run()) { storeModel(p_currentBestModel); break;}
#ifdef STK_MIXTURE_VERBOSE
            stk_cout << "In FullStrategy::run(), Long run Failed.\n";
#endif
    } // end iTry
    delete p_currentBestModel;
    delete p_currentModel;
  } catch (Exception const& e)
  {
    msg_error_ = e.error();
    return false;
  }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << "XemStrategy::run() terminated.\n";
  stk_cout << "-------------------------------\n";
#endif
  if (p_model_->lnLikelihood() <= value)
  {
    msg_error_ += STKERROR_NO_ARG(In XemStrategy::run,No gain\n);
    return false;
  }
  return true;
}


} // namespace STK



