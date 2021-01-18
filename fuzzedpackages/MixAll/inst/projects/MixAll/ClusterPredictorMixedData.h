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
 * created on: 17 Mars 2018construct properly a mixture model.
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file ClusterPredictorMixedData.h
 *  @brief In this file we define the ClusterPredictorMixedData class which
 *  allow to predict class membership for new values in mixed data models
 **/


#ifndef STK_CLUSTERPREDICTORMIXEDDATA_H
#define STK_CLUSTERPREDICTORMIXEDDATA_H

#include "IClusterPredictor.h"

namespace STK
{

/** ClusterPredictorMixedData class allows **/
class ClusterPredictorMixedData: public IClusterPredictor
{
  public:
    /** constructor.
     *  @param model a reference on the current model
     **/
    ClusterPredictorMixedData( Rcpp::S4 model, Rcpp::S4 s4_clusterPredict);
    /** destructor. */
    ~ClusterPredictorMixedData();
    /** run the estimation */
    bool run();

  protected:
    /** list of components of the model from the R side */
    Rcpp::List lcomponent_;
    /** list of data sets of the predictor from the R side */
    Rcpp::List ldata_;
};

} // namespace STK

#endif /* STK_CLUSTERPREDICTORMIXEDDATA_H */
