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
 * Purpose: constants and utilities functions for the project Model.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Model_Util.h
 *  @brief In this file we define the constant and utilities methods used in
 *  the project Model.
 **/

#ifndef STK_MODEL_UTIL_H
#define STK_MODEL_UTIL_H

#include <STKernel/include/STK_String.h>

namespace STK
{

namespace hidden
{
template<class Data> struct ModelTraits;
} // hidden

namespace Model
{

/** @ingroup StatModels
 *  list of the statistical models that can be used
 **/
enum StatisticalModels
{
  ModelBernoulli_pj_,
  ModelDiagGaussian_muj_sj_,
  ModelGamma_aj_bj_,
  MixtureGamma_ajk_bjk_,
  MixtureGamma_ajk_bk_,
  MixtureGamma_ajk_bj_,
  MixtureGamma_ajk_b_,
  MixtureGamma_ak_bjk_,
  MixtureGamma_ak_bk_,
  MixtureGamma_ak_bj_,
  MixtureGamma_ak_b_,
  MixtureGamma_aj_bjk_,
  MixtureGamma_aj_bk_,
  MixtureGamma_a_bjk_,
  MixtureGamma_a_bk_,
  MixtureGaussian_sjk_,
  MixtureGaussian_sk_,
  MixtureGaussian_sj_,
  MixtureGaussian_s_,
  MixtureCategorical_pjk_,
  MixtureCategorical_pk_,
  MixturePoisson_ljk_,
  MixturePoisson_lk_,
  MixturePoisson_ljlk_,
  MixtureKmm_sk_,
  MixtureKmm_s_,
  unknown_model_
};

/** @ingroup StatModels
 *  @brief Model selection criteria
 **/
enum TypeCriterion
{
  unknownCriterion_ = -1, ///< unknown criterion
  aic_,     ///< akaike criterion
  bic_,     ///< information criterion
  icl_,      ///< integrated likelihood criterion
  ml_       ///< Maximum Likelihood criterion
};

/** @ingroup StatModels
 *  convert a String to a TypeReduction.
 *  @param type the String we want to convert
 *  @return the TypeCriterion represented by the String @c type. if the string
 *  does not match any known name, the @c unknown_ type is returned.
 **/
inline TypeCriterion stringToTypeCriterion(const String & type)
{
  if (toUpperString(type) == toUpperString(_T("AIC"))) return aic_;
  if (toUpperString(type) == toUpperString(_T("BIC"))) return bic_;
  if (toUpperString(type) == toUpperString(_T("ICL"))) return icl_;
  if (toUpperString(type) == toUpperString(_T("ML"))) return ml_;
  return unknownCriterion_;
}

/** @ingroup StatModels
 *  convert a TypeCriterion to a String.
 *  @param type the type of criterion we want to convert
 *  @return the string associated to this type.
 **/
inline String typeCriterionToString(const TypeCriterion & type)
{
  switch (type)
  {
    case aic_:
      return String(_T("AIC"));
    case bic_:
      return String(_T("BIC"));
    case icl_:
      return String(_T("ICL"));
    case ml_:
      return String(_T("ML"));
    case unknownCriterion_:
      return String(_T("unknown"));
  }
  return String(_T("unknown"));
}

} // namespace Model

} // namespace STK

#endif /* STK_MODEL_UTIL_H */
