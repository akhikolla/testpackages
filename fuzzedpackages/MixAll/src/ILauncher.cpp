/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

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
 * Project:  MixAll
 * created on: 15 may 2016
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_ILauncher.cpp
 *  @brief In this file we implement the ILauncher which
 *  construct properly a mixture model.
 **/


#include "../inst/projects/MixAll/ILauncher.h"

using namespace Rcpp;

namespace STK
{
/* facade design pattern.
 * The ILauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
ILauncher::ILauncher( Rcpp::S4 model, Rcpp::CharacterVector models)
                    : ILauncherBase(model)
                    , v_models_(models)
{}
/* facade design pattern.
 * The ILauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
ILauncher::ILauncher( Rcpp::S4 model)
                    : ILauncherBase(model)
                    , v_models_()
{}

/* destructor. */
ILauncher::~ILauncher()
{}


/* create the managers for models with real data */
void ILauncher::createContinuousDataSets( String const& idData
                                        , Rcpp::S4 s4_component
                                        , Clust::Mixture model
                                        )
{
  NumericMatrix m_data = s4_component.slot("data");
  handler_.addData(m_data, idData, Clust::mixtureToString(model));
}

/* create the managers for models with real data */
void ILauncher::createDiscreteDataSets( String const& idData
                                      , Rcpp::S4 s4_component
                                      , Clust::Mixture model
                                      )
{
  IntegerMatrix m_data = s4_component.slot("data");
  handler_.addData(m_data, idData, Clust::mixtureToString(model));
}

/* create the mixtures in the given learner */
void ILauncher::createMixtures(IMixtureStatModel* p_model)
{
  p_model->createMixture(diagGaussianManager_);
  p_model->createMixture(poissonManager_);
  p_model->createMixture(gammaManager_);
  p_model->createMixture(categoricalManager_);
}

} // namespace STK

