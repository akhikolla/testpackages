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

/** @file STK_SimpleStrategy.cpp
 *  @brief In this file we implement the SimpleStratefy class
 **/

#include <Sdk.h>

#include "../include/MixtureStrategy/STK_SimpleStrategy.h"
#include "../include/MixtureInit/STK_MixtureInit.h"
#include "../include/MixtureAlgo/STK_MixtureAlgo.h"
#include "../include/STK_IMixtureComposer.h"

namespace STK
{

/* destructor */
SimpleStrategyParam::~SimpleStrategyParam()
{ if (p_algo_) delete p_algo_;}

/* run the simple strategy */
bool SimpleStrategy::run()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("-----------------------------------------------\n");
  stk_cout << _T("Entering SimpleStrategy::run() with: ")
           << _T("nbTry_ = ") << nbTry_ << _T("\n");
#endif
  IMixtureComposer* p_currentModel = 0;
  if (p_model_->state() < 1) { p_model_->randomFuzzyInit();}
  Real value = p_model_->lnLikelihood();
  try
  {
    p_currentModel = p_model_->create();

    // initialize and run algo. break if success
    for (int iTry = 0; iTry < nbTry_; ++iTry)
    {
      // initialize current model
      p_init_->setModel(p_currentModel);
      if (!p_init_->run())
      {
#ifdef STK_MIXTURE_VERBOSE
        stk_cout << "SimpleStrategy::run(), Initialization failed.\n";
#endif
        msg_error_ += STKERROR_NO_ARG(SimpleStrategy::run,Initialization failed\n);
        msg_error_ += p_init_->error();
      } // init step
      p_param_->p_algo_->setModel(p_currentModel);
      // initial
      if (!p_param_->p_algo_->run())
      {
#ifdef STK_MIXTURE_VERBOSE
        stk_cout << "SimpleStrategy::run(), Long run failed.\n";
#endif
        msg_error_ += STKERROR_NO_ARG(SimpleStrategy::run,long run failed\n);
        msg_error_ += p_param_->p_algo_->error();
      }
      else { break;}  // we get a result
    } // iTry
  }
  catch (Exception const& e)
  {
    msg_error_ = e.error();
    return false;
  }
  storeModel(p_currentModel);
  delete p_currentModel;
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << "SimpleStrategy::run() terminated.              \n";
  stk_cout << "-----------------------------------------------\n";
#endif
  if (p_model_->lnLikelihood() <= value)
  {
    msg_error_ += STKERROR_NO_ARG(SimpleStrategy::run,No gain\n);
    return false;
  }
  return true;
}


} // namespace STK



