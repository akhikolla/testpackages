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

#include "../include/MixtureAlgo/STK_MixtureAlgo.h"
#include "../include/STK_IMixtureComposer.h"

namespace STK
{
/* run the CEM algorithm */
bool CEMAlgo::run()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Entering CEMAlgo::run() with:\n")
           << _T("nbIterMax_ = ") << nbIterMax_
           << _T(", epsilon_ = ") << epsilon_ << _T("\n");
#endif
  try
  {
    Real currentLnLikelihood =  p_model_->lnLikelihood();
    int iter;
    for (iter = 0; iter < nbIterMax_; iter++)
    {
      if (p_model_->cStep()<threshold_)
      {
        msg_error_ = STKERROR_NO_ARG(CEMAlgo::run,No more individuals after cStep\n);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("An error occur in CEMAlgo::run():\n") << msg_error_ << _T("\n");
#endif
        return false;
      }
      p_model_->imputationStep();
      p_model_->pStep();
      p_model_->paramUpdateStep();
      Real nb = p_model_->eStep();
      if (nb<threshold_)
      {
        msg_error_ = STKERROR_1ARG(CEMAlgo::run,nb,Not enough individuals after eStep\n);
#ifdef STK_MIXTURE_VERBOSE
        stk_cout << _T("An error occur in CEMAlgo::run():\n") << msg_error_ << _T("\n");
#endif
        return false;
      }
      Real lnLikelihood = p_model_->lnLikelihood();
      if (std::abs(lnLikelihood - currentLnLikelihood) < epsilon_)
      {
#ifdef STK_MIXTURE_VERY_VERBOSE
        stk_cout << _T("Terminating CEMAlgo::run() with:\n")
                 << _T("iter = ") << iter << _T("\n")
                 << _T("delta = ") << lnLikelihood - currentLnLikelihood << _T("\n");
#endif
        break;
      }
      currentLnLikelihood = lnLikelihood;
    }
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("In EMAlgo::run() iteration ") << iter << _T("terminated.\n")
             << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  }
  catch (Clust::exceptions const& error)
  {
    msg_error_ = Clust::exceptionToString(error);
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("An error occur in CEMAlgo::run():\n") << msg_error_ << _T("\n");
#endif
    return false;
  }
  return true;
}

bool EMAlgo::run()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("----------------------------\n");
  stk_cout << _T("Entering EMAlgo::run() with:\n")
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
      p_model_->pStep();
      p_model_->paramUpdateStep();
      Real nb = p_model_->eStep();
      if (nb<threshold_)
      {
        msg_error_ = STKERROR_1ARG(EMAlgo::run,nb,Not enough individuals after eStep\n);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("An error occur in EMAlgo::run():\n") << msg_error_ << _T("\n");
#endif
        return false;
      }
      Real lnLikelihood = p_model_->lnLikelihood();
      // no abs as the likelihood should increase
      if ( (lnLikelihood - currentLnLikelihood) < epsilon_)
      {
#ifdef STK_MIXTURE_VERY_VERBOSE
        stk_cout << _T("Terminating EMAlgo::run() with:\n")
                 << _T("iter = ") << iter << _T("\n")
                 << _T("delta = ") << lnLikelihood - currentLnLikelihood << _T("\n");
#endif
        break;
      }
      currentLnLikelihood = lnLikelihood;
    }
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("In EMAlgo::run() iteration ") << iter << _T("terminated.\n")
             << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  }
  catch (Clust::exceptions const& error)
  {
    msg_error_ = Clust::exceptionToString(error);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("An error occur in EMAlgo::run():\n") << msg_error_ << _T("\n");
#endif
    return false;
  }
  return true;
}

bool SEMAlgo::run()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("-----------------------------\n");
  stk_cout << _T("Entering SEMAlgo::run() with:\n")
           << _T("nbIterMax_ = ") << nbIterMax_ << _T("\n")
           << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  bool result = true;
  try
  {
    int iter;
    for (iter = 0; iter < nbIterMax_; ++iter)
    {
      Real nb = p_model_->sStep(); // simulate labels
      if (nb<threshold_)
      {
        msg_error_ = STKERROR_1ARG(SEMAlgo::run,nb,Not enough individuals after sStep\n);
#ifdef STK_MIXTURE_VERBOSE
        stk_cout << _T("An error occur in SEMAlgo::run():\n") << msg_error_ << _T("\n");
#endif
        result = false;
        break;
      }
      p_model_->samplingStep(); // simulate missing values
      p_model_->pStep();        // estimate proportions
      p_model_->paramUpdateStep();        // estimate parameters
      nb = p_model_->eStep();   // update tik and lnLikelihood
      if (nb<threshold_)
      {
        msg_error_ = STKERROR_1ARG(SEMAlgo::run,nb,Not enough individuals after eStep\n);
#ifdef STK_MIXTURE_VERBOSE
        stk_cout << _T("An error occur in SEMAlgo::run():\n") << msg_error_ << _T("\n");
#endif
        return false;
        break;
      }
      p_model_->storeIntermediateResults(iter+1); // store current parameters
    }
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("In SEMAlgo::run() iterations terminated.\n")
             << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  }
  catch (Clust::exceptions const& error)
  {
    msg_error_ = Clust::exceptionToString(error);
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("An error occur in SEMAlgo::run(): ") << msg_error_ << _T("\n");
#endif
    result = false;
  }
  if (result)
  {
    // set averaged parameters
    p_model_->setParametersStep();
#ifdef STK_MIXTURE_VERY_VERBOSE
    stk_cout << _T("\nIn SEMAlgo::run(), setParameters done.\n")
             << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  }
  else
  { p_model_->releaseIntermediateResults();}
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Terminating SEMAlgo::run()\n");
  stk_cout << _T("--------------------------\n");
#endif
  return result;
}

bool SemiSEMAlgo::run()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("---------------------------------\n");
  stk_cout << _T("Entering SemiSEMAlgo::run() with:\n")
           << _T("nbIterMax_ = ") << nbIterMax_ << _T("\n")
           << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  bool result = true;
  try
  {
    int iter;
    for (iter = 0; iter < nbIterMax_; ++iter)
    {
      p_model_->samplingStep();
      p_model_->pStep();
      p_model_->paramUpdateStep();
      Real nb = p_model_->eStep();
      if (nb<threshold_)
      {
        msg_error_ = STKERROR_1ARG(SemiSEMAlgo::run,nb,Not enough individuals after eStep\n);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("An exception occur in SemiSEMAlgo::run(): ") << msg_error_ << _T("\n");
#endif
        result = false;
        break;
      }
      p_model_->storeIntermediateResults(iter+1); // store current parameters
    }
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("In SemiSEMAlgo::run() iteration ") << iter << _T("terminated.\n")
             << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  }
  catch (Clust::exceptions const& error)
  {
    msg_error_ = Clust::exceptionToString(error);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("An exception occur in SemiSEMAlgo::run(): ") << msg_error_ << _T("\n");
#endif
    result = false;
  }
  // set averaged parameters
  if (result)
  {
    p_model_->setParametersStep();
#ifdef STK_MIXTURE_VERY_VERBOSE
    stk_cout << _T("\nIn SemiSEMAlgo::run(), setParameters done.\n")
             << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  }
  else // no result
  { p_model_->releaseIntermediateResults();}
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Terminating SemiSEMAlgo::run()\n");
  stk_cout << _T("------------------------------\n");
#endif
  return result;
}


} // namespace STK
