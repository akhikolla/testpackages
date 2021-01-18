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

/** @file ClusterPredictor.h
 *  @brief In this file we define the ClusterPredictor class which
 *  allow to predict class membership for new values
 **/


#ifndef STK_CLUSTERPREDICTOR_H
#define STK_CLUSTERPREDICTOR_H

#include "IClusterPredictor.h"

namespace STK
{

/** ClusterPredictor class allows **/
class ClusterPredictor: public IClusterPredictor
{
  public:
    /** constructor.
     *  @param model a reference on the current model
     **/
    ClusterPredictor( Rcpp::S4 model, Rcpp::S4 s4_clusterPredict);
    /** destructor. */
    ~ClusterPredictor();
    /** run the estimation */
    bool run();

  protected:
    /** component of the model from the R side */
    Rcpp::S4 s4_component_;
};

} // namespace STK

#endif /* STK_CLUSTERPREDICTOR_H */
