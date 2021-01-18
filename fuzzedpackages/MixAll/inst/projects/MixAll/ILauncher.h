/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2015  Serge Iovleff, University Lille 1, Inria

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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 10 May 2016
 * Author:   Iovleff, serge.iovleff@stkpp.org
 **/

/** @file ILauncher.h
 *  @brief In this file we define the interface class for launchers.
 **/


#ifndef STK_ILAUNCHER_H
#define STK_ILAUNCHER_H

#include "ILauncherBase.h"

namespace STK
{

/** The ILauncher allow to create the composer or learner for estimate or
 *  learn a  mixture model with less effort.
 **/
class ILauncher: public ILauncherBase
{
  public:
    /** constructor.
     * @param model a reference on the current model
     * @param models a list of model name
     **/
    ILauncher( Rcpp::S4 model, Rcpp::CharacterVector models);
    /** constructor with a list of component.
     *  @param model a reference on the current model
     **/
    ILauncher( Rcpp::S4 model);
    /** destructor. */
    virtual ~ILauncher();
    /** @return the model */
    inline Rcpp::S4 const& s4_model() const { return s4_model_;}

  protected:
    /** create data sets */
    template<int Rtype>
    void createDataSets(Rcpp::Matrix<Rtype> const& data, String const& idData, Clust::Mixture model)
    { handler_.addData(data, idData, Clust::mixtureToString(model));}

    /** create the data sets with real data */
    void createContinuousDataSets(String const& idData, Rcpp::S4 s4_component, Clust::Mixture model);
    /** create the data sets with integer data */
    void createDiscreteDataSets(String const& idData, Rcpp::S4 s4_component, Clust::Mixture model);
    /** create the mixtures in the given model */
    void createMixtures(IMixtureStatModel* p_model);

    /** vector with the model names to try */
    Rcpp::CharacterVector v_models_;
};

} // namespace STK

#endif /* STK_ILAUNCHER_H */
