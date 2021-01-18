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
 * Project:  stkpp::Model
 * created on: 22 juil. 2011
 * Purpose: define the Interface base class IMixtureCriterion.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_IMixtureCriterion.h
 *  @brief In this file we define the interface bas class for computing penalized
 *  criterion on mixture models
 **/

#ifndef STK_IMIXTURECRITERION_H
#define STK_IMIXTURECRITERION_H

#include <Sdk.h>
#include "../STK_IMixtureStatModel.h"

namespace STK
{

/** @ingroup Clustering
 * @brief Interface base class for the selection model criterion. The pure
 * virtual function @c run will be implemented in derived class and compute
 * the value_ member.
  */
class IMixtureCriterion: public IRunnerBase
{
  protected:
    /** Default Constructor. */
    inline IMixtureCriterion() : p_composer_(), value_(Arithmetic<Real>::NA()) {}
    /** Constructor.
     *  @param p_composer a pointer on the current model
     **/
    inline IMixtureCriterion( IMixtureStatModel const* p_composer)
                            : p_composer_(p_composer), value_(Arithmetic<Real>::NA()){}
    /** copy Constructor.
     *  @param criterion the criterion to copy
     **/
    inline IMixtureCriterion( IMixtureCriterion const& criterion)
                            : p_composer_(criterion.p_composer_)
                            , value_(criterion.value_) {}
  public:
    /** Destructor */
    inline virtual ~IMixtureCriterion() {}
    /** @return The value of the criterion */
    inline Real const& value() const { return value_;}
    /** @param p_composer a pointer on the current model to set */
    inline void setModel( IMixtureStatModel const* p_composer)
    { p_composer_ = p_composer;}

  protected:
    /** The current statistical model to use*/
    IMixtureStatModel const* p_composer_;
    /** Computed value of the criterion */
    Real value_;
};


} // namespace STK

#endif /** STK_MIXTURECRITERION_H */
