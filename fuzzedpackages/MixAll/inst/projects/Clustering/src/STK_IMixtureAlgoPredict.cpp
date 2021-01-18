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

/** @file STK_IMixtureAlgoPredict.cpp
 *  @brief In this file we implement the run method of the interface prediction algorithms.
 **/

#include "../include/MixtureAlgo/STK_IMixtureAlgoPredict.h"
#include "../include/STK_IMixtureComposer.h"

namespace STK
{
/* default constructor */
IMixtureAlgoPredict::IMixtureAlgoPredict(): IRunnerBase()
                                          , p_model_(0)
                                          , nbIterBurn_(0), nbIterLong_(0), epsilon_(0.) {}

/* Copy constructor.
 *  @param algo the algorithm to copy */
IMixtureAlgoPredict::IMixtureAlgoPredict( IMixtureAlgoPredict const& algo): IRunnerBase(algo)
                                         , p_model_(algo.p_model_)
                                         , nbIterBurn_(algo.nbIterBurn_)
                                         , nbIterLong_(algo.nbIterLong_)
                                         , epsilon_(algo.epsilon_)
{}

/* destructor */
IMixtureAlgoPredict::~IMixtureAlgoPredict() {}

// threshold_ is set to this value in order to get stability on the results
void IMixtureAlgoPredict::setModel(IMixtureComposer* p_model)
{ p_model_ = p_model;}


/* predict class labels when there is no missing values
 * @return @c true if no error occur, @c false otherwise
 */
bool IMixtureAlgoPredict::predictBayesClassifier()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("-------------------------------\n");
  stk_cout << _T("Entering IMixtureAlgoPredict::predictBayesClassifier()\n");
#endif


  try
  {
    p_model_->initializeStep();
    p_model_->eStep();
    p_model_->mapStep();
    p_model_->finalizeStep();
  }
  catch (Exception const& error)
  {
    msg_error_ = error.error();
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("An error occur in IMixtureAlgoPredict::run():\n") << msg_error_ << _T("\n");
#endif
    return false;
  }
  catch (Clust::exceptions const& error)
  {
    msg_error_ = Clust::exceptionToString(error);
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("An error occur in IMixtureAlgoPredict::run():\n") << msg_error_ << _T("\n");
#endif
    return false;
  }
  return true;
}

// burn step
bool IMixtureAlgoPredict::burnStep()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("-------------------------------\n");
  stk_cout << _T("Entering IMixtureAlgoPredict::burnStep()\n");
#endif
  try
  {
    p_model_->initializeStep();
    for (int iter = 0; iter < nbIterLong_; ++iter)
    {
      p_model_->sStep(); // simulate labels
      p_model_->samplingStep(); // simulate missing values
      p_model_->eStep();   // update tik and lnLikelihood
      p_model_->storeIntermediateResults(iter+1); // store current parameters
    }
  }
  catch (Clust::exceptions const& error)
  {
    msg_error_ = Clust::exceptionToString(error);
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("An error occur in IMixtureAlgoPredict::run(): ") << msg_error_ << _T("\n");
#endif
    p_model_->setParametersStep(); // compute current parameters value
    return false;
  }
  p_model_->setParametersStep(); // compute current parameters value
  return true;
}

} // namespace STK
