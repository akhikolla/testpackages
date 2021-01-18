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

/** @file STK_MixtureCriterion.h
 *  @brief In this file we define the classes for computing the penalized
 *  criterion on mixture models
 **/

#ifndef STK_MIXTURECRITERION_H
#define STK_MIXTURECRITERION_H

#include "STK_IMixtureCriterion.h"

namespace STK
{

/** @ingroup Clustering
 *  @brief Derived class of Criterion for computing the AIC Criterion
 *  The AIC criteria of a given model M is a penalization of the likelihood
 *  given by the formula
 *  \f[
 *  -2 \cdot \ln{p(x|M)} \approx \mathrm{AIC} = {-2 \cdot \ln{L} +2 D }
 *  \f]
 *  where \f$ L \f$ represents the likelihood of the observations and \f$ D \f$
 *  the number of free parameter of the model.
 **/
class AICMixtureCriterion: public IMixtureCriterion
{
  public:
    /** Default Constructor. */
    inline AICMixtureCriterion(): IMixtureCriterion() {}
    /** Constructor.
     *  @param p_composer apointer on the current model
     **/
    inline AICMixtureCriterion( IMixtureStatModel* const p_composer)
                              : IMixtureCriterion(p_composer) {}
    /** copy Constructor.
     *  @param criterion the criterion to copy
     **/
    inline AICMixtureCriterion( AICMixtureCriterion const& criterion)
                              : IMixtureCriterion(criterion) {}
    /** virtual destructor. */
    inline virtual ~AICMixtureCriterion() {}
    /** clone pattern */
    inline AICMixtureCriterion* clone() const { return new AICMixtureCriterion(*this);}
    /** implementation of the virtual method run */
    virtual bool run();
};


/** @ingroup Clustering
 *  @brief Derived class of Criterion for computing the BIC Criterion
 *  The Bic criteria of a model M is a penalization of the likelihood given by
 *  the formula
 *  \f[
 *  -2 \cdot \ln{p(x|M)} \approx \mathrm{BIC} = {-2 \cdot \ln{L} + D \ln(n) }
 *  \f]
 *  where \f$ L \f$ represents the likelihood of the observations, \f$ D \f$ the
 *  number of free parameter of the model and \f$ n \f$ the number of sample.
 **/
class BICMixtureCriterion: public IMixtureCriterion
{
  public:
    /** Default Constructor. */
    inline BICMixtureCriterion(): IMixtureCriterion() {}
    /** Constructor.
     *  @param p_composer apointer on the current model
     **/
    inline BICMixtureCriterion( IMixtureStatModel* const p_composer)
                              : IMixtureCriterion(p_composer) {}
    /** copy Constructor.
     *  @param criterion the criterion to copy
     **/
    inline BICMixtureCriterion( BICMixtureCriterion const& criterion)
                              : IMixtureCriterion(criterion) {}
    /** virtual destructor. */
    inline virtual ~BICMixtureCriterion() {}
    /** clone pattern */
    inline BICMixtureCriterion* clone() const { return new BICMixtureCriterion(*this);}
    /** implementation of the virtual method run */
    virtual bool run();
};

/** @ingroup Clustering
 *  @brief Derived class of IMixtureCriterion for computing the ICL Criterion
 *  The ICL criteria of a model M is a penalization of the likelihood given by the formula
 *  \f[
 *  -2 \ln{p(x|M)} \approx \mathrm{ICL} = { -2 \ln{L} + D \log(n) - 2 \sum_i\sum_k t_{ik}\log(t_{ik}) }
 *  \f]
 *  where \f$ L \f$ represents the likelihood of the observations and \f$ D \f$
 *  the number of free parameters of the model.
 **/
class ICLMixtureCriterion: public IMixtureCriterion
{
  public:
    /** Default Constructor. */
    inline ICLMixtureCriterion(): IMixtureCriterion() {}
    /** Constructor.
     *  @param p_composer a pointer on the current composer
     **/
    inline ICLMixtureCriterion( IMixtureStatModel const* p_composer)
                              : IMixtureCriterion(p_composer) {}
    /** copy Constructor.
     *  @param criterion the criterion to copy
     **/
    inline ICLMixtureCriterion( ICLMixtureCriterion const& criterion)
                              : IMixtureCriterion(criterion) {}
    /** virtual destructor. */
    inline virtual ~ICLMixtureCriterion() {}
    /** clone pattern */
    inline ICLMixtureCriterion* clone() const { return new ICLMixtureCriterion(*this);}
    /** implementation of the virtual method run */
    virtual bool run();
};

/** @ingroup Clustering
 *  @brief Derived class of Criterion for computing the Maximum Likelihood Criterion
 *  This criterion does not penalize the model and return minux two times the
 *  log likelihood of the model.
 *  \f[
 *    -2 \ln{L}
 *  \f]
 *
 **/
class MLMixtureCriterion: public IMixtureCriterion
{
  public:
    /** Default Constructor. */
    inline MLMixtureCriterion(): IMixtureCriterion() {}
    /** Constructor.
     *  @param p_composer apointer on the current model
     **/
    inline MLMixtureCriterion( IMixtureStatModel* const p_composer)
                              : IMixtureCriterion(p_composer) {}
    /** copy Constructor.
     *  @param criterion the criterion to copy
     **/
    inline MLMixtureCriterion( MLMixtureCriterion const& criterion)
                              : IMixtureCriterion(criterion) {}
    /** virtual destructor. */
    inline virtual ~MLMixtureCriterion() {}
    /** clone pattern */
    inline MLMixtureCriterion* clone() const { return new MLMixtureCriterion(*this);}
    /** implementation of the virtual method run */
   virtual bool run();
};


} // namespace STK

#endif /** STK_MIXTURECRITERION_H */
