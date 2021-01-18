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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 14 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/


/**@file STK_IMixture.h
 * @brief define the main interface for linking specific mixture model to the
 * composer.
 */

#include "../include/STK_IMixture.h"
#include "../include/STK_IMixtureStatModel.h"

namespace STK
{

/* constructor */
IMixture::IMixture( String const& idData)
                  : p_composer_(0), idData_(idData)
{}

/* copy constructor */
IMixture::IMixture( IMixture const& mixture)
                  : p_composer_(0)
                  , idData_(mixture.idData_)
{}
/* Virtual destructor. */
IMixture::~IMixture(){}

/* set the mixture composer to the mixture */
void IMixture::setMixtureModel( IMixtureStatModel const* p_composer) { p_composer_ = p_composer;}

/* This function can be used in derived classes to get number of samples.
 *  @return Number of samples.
 *  @note STK_IMixtureStatModel.h is include only here
 */
int IMixture::nbSample() const { return p_composer_->nbSample();}
/* This function can be used in derived classes to get number of classes.
 *  @return Number of classes.
 */
int IMixture::nbCluster() const { return p_composer_->nbCluster();}

/* This function can be used in derived classes to get proportions from the framework.
 * @return Pointer to proportions.
 */
CPointX const* IMixture::p_pk() const { return &(p_composer_->pk());}
/* This function can be used in derived classes to get proportions from the framework.
 * @return Pointer to proportions.
 */
CPointX const* IMixture::p_tk() const { return &(p_composer_->tk());}
/* This function can be used in derived classes to get proportions from the framework.
 * @return Pointer to proportions.
 */
CArrayXX const* IMixture::p_tik() const { return &(p_composer_->tik());}
/* This function can be used in derived classes to get proportions from the framework.
 * @return Pointer to proportions.
 */
CVectorXi const* IMixture::p_zi() const { return &(p_composer_->zi());}

} // namespace STK

