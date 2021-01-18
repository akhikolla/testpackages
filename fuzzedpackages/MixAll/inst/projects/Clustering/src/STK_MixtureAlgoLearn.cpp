/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_MixtureAlgo.cpp
 *  @brief In this file we implement the run method of the mixture algorithms.
 **/

#include <Sdk.h>

#include "../include/MixtureAlgo/STK_MixtureAlgoLearn.h"
#include "../include/STK_IMixtureLearner.h"

namespace STK
{

bool ImputeAlgo::run()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("--------------------------------\n");
  stk_cout << _T("Entering ImputeAlgo::run() with:\n")
           << _T("nbIterMax_ = ") << nbIterMax_ << _T("\n")
           << _T("epsilon_ = ") << epsilon_ << _T("\n");
#endif

  try
  {
    Real currentLnLikelihood = p_model_->lnLikelihood();
    int iter;
    for (iter = 0; iter < nbIterMax_; iter++)
    {
      p_model_->imputationStep();
      p_model_->paramUpdateStep();
      Real lnLikelihood = p_model_->lnLikelihood();
      // no abs as the likelihood should increase
      if ( (lnLikelihood - currentLnLikelihood) < epsilon_)
      {
#ifdef STK_MIXTURE_VERY_VERBOSE
        stk_cout << _T("Terminating ImputeAlgo::run() with:\n")
                 << _T("iter = ") << iter << _T("\n")
                 << _T("delta = ") << lnLikelihood - currentLnLikelihood << _T("\n");
#endif
        break;
      }
      currentLnLikelihood = lnLikelihood;
    }
    // finalize and compute posterior probabilities and predicted values
    p_model_->finalizeStep();
    p_model_->mapStep();
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("In ImputeAlgo::run() iteration ") << iter << _T("terminated.\n")
             << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  }
  catch (Clust::exceptions const& error)
  {
    msg_error_ = Clust::exceptionToString(error);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("An error occur in ImputeAlgo::run():\n") << msg_error_ << _T("\n");
#endif
    return false;
  }
  return true;
}

bool SimulAlgo::run()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("-------------------------------\n");
  stk_cout << _T("Entering SimulAlgo::run() with:\n")
           << _T("nbIterMax_ = ") << nbIterMax_ << _T("\n")
           << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Parameters of the model\n");
  p_model_->writeParameters(stk_cout);
#endif
  bool result = true;
  try
  {
    int iter;
    for (iter = 0; iter < nbIterMax_; ++iter)
    {
      p_model_->samplingStep(); // simulate missing values
      p_model_->paramUpdateStep(); // estimate parameters
      p_model_->storeIntermediateResults(iter+1); // store current parameters
    }
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("In SimulAlgo::run() iterations terminated.\n")
             << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Parameters of the model\n");
  p_model_->writeParameters(stk_cout);
#endif
  }
  catch (Clust::exceptions const& error)
  {
    msg_error_ = Clust::exceptionToString(error);
    p_model_->finalizeStep();
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("An error occur in SimulAlgo::run(): ") << msg_error_ << _T("\n");
#endif
    result = false;
  }
  if (result)
  {
    // set averaged parameters
    p_model_->finalizeStep();
    p_model_->mapStep();
#ifdef STK_MIXTURE_VERY_VERBOSE
    stk_cout << _T("\nIn SimulAlgo::run(), setParameters done.\n")
             << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  }
  else
  { p_model_->releaseIntermediateResults();}
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Terminating SimulAlgo::run()\n");
  stk_cout << _T("----------------------------\n");
#endif
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Parameters of the model\n");
  p_model_->writeParameters(stk_cout);
#endif
  return result;
}

} // namespace STK
