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

/** @file STK_IMixtureStrategy.cpp
 *  @brief In this file we implement the interface base class for estimation strategies
 **/

#include "../include/STK_IMixtureComposer.h"
#include "../include/MixtureInit/STK_MixtureInit.h"
#include "../include/MixtureStrategy/STK_IMixtureStrategy.h"

namespace STK
{

/* default constructor.
 *  @param p_model the model to estimate
 **/
IMixtureStrategy::IMixtureStrategy( IMixtureComposer*& p_model)
                                  : IRunnerBase(), nbTry_(1), p_model_(p_model), p_init_(0)
{}
/* copy constructor
 *  @param strategy the strategy to copy
 **/
IMixtureStrategy::IMixtureStrategy( IMixtureStrategy const& strategy)
                                  : IRunnerBase(strategy), nbTry_(strategy.nbTry_)
                                  , p_model_(strategy.p_model_)
                                  , p_init_(strategy.p_init_->clone())
{}

/* store a model in p_model_ if it is better.
 * @param p_otherModel the model to store
 **/
void IMixtureStrategy::storeModel(IMixtureComposer*& p_otherModel)
{ if (p_model_->lnLikelihood()<p_otherModel->lnLikelihood())
  { std::swap(p_model_, p_otherModel);}
}

/* destructor */
IMixtureStrategy::~IMixtureStrategy() { if (p_init_) delete p_init_;}


} // namespace STK



