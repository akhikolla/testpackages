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
 * created on: 4 sept. 2013
 * Author:   iovleff, s...DOTi...ATstkppDOTorg
 **/

/** @file STK_Clust_Util.cpp
 *  @brief In this file we implement the utilities functions of the Clustering project.
 **/
#include "../include/MixtureCriterion/STK_MixtureCriterion.h"
#include "../include/MixtureInit/STK_MixtureInit.h"
#include "../include/MixtureStrategy/STK_MixtureStrategy.h"
#include "../include/STK_IMixtureComposer.h"
#include "../include/MixtureAlgo/STK_MixtureAlgo.h"
#include "../include/MixtureAlgo/STK_MixtureAlgoLearn.h"
#include "../include/MixtureAlgo/STK_MixtureAlgoPredict.h"

namespace STK
{

namespace Clust
{

/* @ingroup Clustering
 *  Convert a String to a initType. The recognized strings are
 * <table>
 * <tr> <th> Initialization   </th></tr>
 * <tr> <td> "randomInit"</td></tr>
 * <tr> <td> "randomParamInit"</td></tr>
 * <tr> <td> "randomClassInit"</td></tr>
 * <tr> <td> "randomFuzzyInit"</td></tr>
 * <tr> <td> "random"         </td></tr>
 * <tr> <td> "class"          </td></tr>
 * <tr> <td> "fuzzy"          </td></tr>
 * </table>
 *  @param type the type of initialization wanted
 *  @return the initType corresponding (default is randomClassInit)
 *  @note if the string is not found in the list above,the type Clust::randomClassInit_
 *  is returned.
 **/
initType stringToInit( String const& type)
{
  if (toUpperString(type) == toUpperString(_T("randomInit"))) return randomParamInit_;
  if (toUpperString(type) == toUpperString(_T("randomParamInit"))) return randomParamInit_;
  if (toUpperString(type) == toUpperString(_T("randomClassInit"))) return randomClassInit_;
  if (toUpperString(type) == toUpperString(_T("randomFuzzyInit"))) return randomFuzzyInit_;
  if (toUpperString(type) == toUpperString(_T("random"))) return randomParamInit_;
  if (toUpperString(type) == toUpperString(_T("class"))) return randomClassInit_;
  if (toUpperString(type) == toUpperString(_T("fuzzy"))) return randomFuzzyInit_;
  return randomClassInit_;
}

/* @ingroup Clustering
 *  Convert a String to an algoType. The recognized strings are
 * <table>
 * <tr> <th> Algorithm   </th></tr>
 * <tr> <td> "emAlgo"</td></tr>
 * <tr> <td> "cemAlgo"</td></tr>
 * <tr> <td> "semAlgo"</td></tr>
 * <tr> <td> "em"</td></tr>
 * <tr> <td> "cem"         </td></tr>
 * <tr> <td> "sem"          </td></tr>
 * </table>
 *  @param type the type of algorithm wanted
 *  @return the algoType corresponding (default is emAlgo)
 *  @note if the string is not found in the list above,the type Clust::emAlgo_
 *  is returned.
 **/
algoType stringToAlgo( String const& type)
{
  if (toUpperString(type) == toUpperString(_T("emAlgo"))) return emAlgo_;
  if (toUpperString(type) == toUpperString(_T("cemAlgo"))) return cemAlgo_;
  if (toUpperString(type) == toUpperString(_T("semAlgo"))) return semAlgo_;
  if (toUpperString(type) == toUpperString(_T("semiSemAlgo"))) return semiSemAlgo_;
  if (toUpperString(type) == toUpperString(_T("em"))) return emAlgo_;
  if (toUpperString(type) == toUpperString(_T("cem"))) return cemAlgo_;
  if (toUpperString(type) == toUpperString(_T("sem"))) return semAlgo_;
  if (toUpperString(type) == toUpperString(_T("semiSem"))) return semiSemAlgo_;
  return emAlgo_;
}

/* @ingroup Clustering
 *  Convert a String to an algoPredictType. The recognized strings are
 * <table>
 * <tr> <th> Algorithm     </th></tr>
 * <tr> <td> "em"          </td></tr>
 * <tr> <td> "semiSem"         </td></tr>
 * </table>
 *  @param type the type of algorithm wanted
 *  @return the algoPredictType corresponding (default is em)
 *  @note The capitalized letters have no effect and if the string is not found
 *  in the list above, the type Clust::emPredictAlgo_ is returned.
 **/
algoPredictType stringToPredictAlgo( String const& type)
{
  if (toUpperString(type) == toUpperString(_T("em"))) return emPredictAlgo_;
  if (toUpperString(type) == toUpperString(_T("semiSem"))) return semiSEMPredictAlgo_;
  return emPredictAlgo_;
}

/* @ingroup Clustering
 *  Convert a String to an algoLearnType. The recognized strings are
 * <table>
 * <tr> <th> Algorithm     </th></tr>
 * <tr> <td> "imputeAlgo"  </td></tr>
 * <tr> <td> "simulAlgo"   </td></tr>
 * <tr> <td> "impute"      </td></tr>
 * <tr> <td> "simul"       </td></tr>
 * </table>
 *  @param type the type of algorithm wanted
 *  @return the algoType corresponding (default is emAlgo)
 *  @note The capitalized letters have no effect and if the string is not found
 *  in the list above,the type Clust::emAlgo_ is returned.
 **/
algoLearnType stringToLearnAlgo( String const& type)
{
  if (toUpperString(type) == toUpperString(_T("imputeAlgo"))) return imputeAlgo_;
  if (toUpperString(type) == toUpperString(_T("simulAlgo"))) return simulAlgo_;
  if (toUpperString(type) == toUpperString(_T("impute"))) return imputeAlgo_;
  if (toUpperString(type) == toUpperString(_T("simul"))) return simulAlgo_;
  return imputeAlgo_;
}



/* @ingroup Clustering
 *  Convert a String to an criterionType. The recognized strings are
 * <table>
 * <tr> <th> Criterion </th></tr>
 * <tr> <td> "AIC"     </td></tr>
 * <tr> <td> "BIC"     </td></tr>
 * <tr> <td> "ICL      </td></tr>
 * <tr> <td> "ML"      </td></tr>
 * </table>
 *  @param type the type of criterion wanted
 *  @return the criterionType corresponding
 *  @note The capitalized letters have no effect and if the string is not found
 *  in the list above,the type Clust::bic_ is returned.
 **/
criterionType stringToCriterion( String const& type)
{
  if (toUpperString(type) == toUpperString(_T("AIC"))) return aic_;
  if (toUpperString(type) == toUpperString(_T("BIC"))) return bic_;
  if (toUpperString(type) == toUpperString(_T("ICL"))) return icl_;
  if (toUpperString(type) == toUpperString(_T("ML"))) return ml_;
  return bic_;
}

/* @ingroup Clustering
 *  convert a TypeReduction to a String.
 *  @param type the type of exception that occur
 *  @return the string associated to this type.
 **/
String exceptionToString( exceptions const& type)
{
  if (type == randomInitFail_)      return String(_T("RandomInit fail"));
  if (type == randomClassInitFail_) return String(_T("RandomClassInit fail"));
  if (type == randomFuzzyInitFail_) return String(_T("RandomFuzzyInit fail"));
  if (type == randomParamInitFail_) return String(_T("RandomParamInit fail"));
  if (type == initializeStepFail_)  return String(_T("initializeStep fail"));
  if (type == estimFail_)   return String(_T("Estimation fail"));
  if (type == mStepFail_)   return String(_T("mStep fail"));
  if (type == eStepFail_)   return String(_T("eStep fail"));
  if (type == mapStepFail_) return String(_T("mapStep fail"));
  if (type == cStepFail_)   return String(_T("cStep fail"));
  if (type == sStepFail_)   return String(_T("sStep fail"));
  return String(_T("unknown exception"));
}

/* @ingroup Clustering
 *  convert a Mixture to a MixtureClass.
 *  @param type the type of Mixture
 *  @return the MixtureClass associated to this Mixture.
 **/
MixtureClass mixtureToMixtureClass( Mixture const& type)
{
  if (type == Gamma_ajk_bjk_)     return Gamma_;
  if (type == Gamma_ajk_bk_)      return Gamma_;
  if (type == Gamma_ajk_bj_)      return Gamma_;
  if (type == Gamma_ajk_b_)       return Gamma_;
  if (type == Gamma_ak_bjk_)      return Gamma_;
  if (type == Gamma_ak_bk_)       return Gamma_;
  if (type == Gamma_ak_bj_)       return Gamma_;
  if (type == Gamma_ak_b_)        return Gamma_;
  if (type == Gamma_aj_bjk_)      return Gamma_;
  if (type == Gamma_aj_bk_)       return Gamma_;
  if (type == Gamma_a_bjk_)       return Gamma_;
  if (type == Gamma_a_bk_)        return Gamma_;
  if (type == Gaussian_sjk_)      return DiagGaussian_;
  if (type == Gaussian_sk_)       return DiagGaussian_;
  if (type == Gaussian_sj_)       return DiagGaussian_;
  if (type == Gaussian_s_)        return DiagGaussian_;
  if (type == Gaussian_sjsk_)     return DiagGaussian_;
  if (type == HDGaussian_AjkBkQkDk_) return HDGaussian_;
  if (type == HDGaussian_AjkBkQkD_)  return HDGaussian_;
  if (type == HDGaussian_AjkBkQDk_)  return HDGaussian_;
  if (type == HDGaussian_AjkBkQD_)   return HDGaussian_;
  if (type == HDGaussian_AjkBQkDk_)  return HDGaussian_;
  if (type == HDGaussian_AjkBQkD_)   return HDGaussian_;
  if (type == HDGaussian_AjkBQDk_)   return HDGaussian_;
  if (type == HDGaussian_AjkBQD_)    return HDGaussian_;
  if (type == HDGaussian_AkBkQkDk_)  return HDGaussian_;
  if (type == HDGaussian_AkBkQkD_)   return HDGaussian_;
  if (type == HDGaussian_AkBkQDk_)   return HDGaussian_;
  if (type == HDGaussian_AkBkQD_)    return HDGaussian_;
  if (type == HDGaussian_AkBQkDk_)   return HDGaussian_;
  if (type == HDGaussian_AkBQkD_)    return HDGaussian_;
  if (type == HDGaussian_AkBQDk_)    return HDGaussian_;
  if (type == HDGaussian_AkBQD_)     return HDGaussian_;
  if (type == HDGaussian_AjBkQkDk_)  return HDGaussian_;
  if (type == HDGaussian_AjBkQkD_)   return HDGaussian_;
  if (type == HDGaussian_AjBkQDk_)   return HDGaussian_;
  if (type == HDGaussian_AjBkQD_)    return HDGaussian_;
  if (type == HDGaussian_AjBQkDk_)   return HDGaussian_;
  if (type == HDGaussian_AjBQkD_)    return HDGaussian_;
  if (type == HDGaussian_ABkQkDk_)   return HDGaussian_;
  if (type == HDGaussian_ABkQkD_)    return HDGaussian_;
  if (type == HDGaussian_ABkQDk_)    return HDGaussian_;
  if (type == HDGaussian_ABkQD_)     return HDGaussian_;
  if (type == HDGaussian_ABQkDk_)    return HDGaussian_;
  if (type == HDGaussian_ABQkD_)     return HDGaussian_;
  if (type == Categorical_pjk_)   return Categorical_;
  if (type == Categorical_pk_)    return Categorical_;
  if (type == Poisson_ljk_)       return Poisson_;
  if (type == Poisson_lk_)        return Poisson_;
  if (type == Poisson_ljlk_)      return Poisson_;
  if (type == Kmm_sk_) return Kmm_;
  if (type == Kmm_s_)  return Kmm_;
  if (type == unknown_mixture_)   return unknown_mixture_class_;
  return unknown_mixture_class_; // avoid compiler warning
}

/* @ingroup Clustering
 *  convert a String to an Mixture.
 *  @param type the String we want to convert
 *  @return the Mixture represented by the String @c type. if the string
 *  does not match any known name, the @c unknown_mixture_ type is returned.
 **/
Mixture stringToMixture( String const& type)
{
  if (toUpperString(type) == toUpperString(_T("Gamma_ajk_bjk")))    return Gamma_ajk_bjk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_ajk_bk")))     return Gamma_ajk_bk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_ajk_bj")))     return Gamma_ajk_bj_;
  if (toUpperString(type) == toUpperString(_T("Gamma_ajk_b")))      return Gamma_ajk_b_;
  if (toUpperString(type) == toUpperString(_T("Gamma_ak_bjk")))     return Gamma_ak_bjk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_ak_bk")))      return Gamma_ak_bk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_ak_bj")))      return Gamma_ak_bj_;
  if (toUpperString(type) == toUpperString(_T("Gamma_ak_b")))       return Gamma_ak_b_;
  if (toUpperString(type) == toUpperString(_T("Gamma_aj_bjk")))     return Gamma_aj_bjk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_aj_bk")))      return Gamma_aj_bk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_a_bjk")))      return Gamma_a_bjk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_a_bk")))       return Gamma_a_bk_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_sjk")))     return Gaussian_sjk_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_sk")))      return Gaussian_sk_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_sj")))      return Gaussian_sj_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_s")))       return Gaussian_s_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_sjsk")))    return Gaussian_sjsk_;
  if (toUpperString(type) == toUpperString(_T("Categorical_pjk")))  return Categorical_pjk_;
  if (toUpperString(type) == toUpperString(_T("Categorical_pk")))   return Categorical_pk_;
  if (toUpperString(type) == toUpperString(_T("Poisson_ljk")))      return Poisson_ljk_;
  if (toUpperString(type) == toUpperString(_T("Poisson_lk")))       return Poisson_lk_;
  if (toUpperString(type) == toUpperString(_T("Poisson_ljlk")))     return Poisson_ljlk_;
  if (toUpperString(type) == toUpperString(_T("Kmm_sk")))return Kmm_sk_;
  if (toUpperString(type) == toUpperString(_T("Kmm_s"))) return Kmm_s_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjkBkQkDk")))return HDGaussian_AjkBkQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjkBkQkD"))) return HDGaussian_AjkBkQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjkBkQDk"))) return HDGaussian_AjkBkQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjkBkQD")))  return HDGaussian_AjkBkQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjkBQkDk"))) return HDGaussian_AjkBQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjkBQkD")))  return HDGaussian_AjkBQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjkBQDk")))  return HDGaussian_AjkBQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjkBQD")))   return HDGaussian_AjkBQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AkBkQkDk"))) return HDGaussian_AkBkQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AkBkQkD")))  return HDGaussian_AkBkQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AkBkQDk")))  return HDGaussian_AkBkQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AkBkQD")))   return HDGaussian_AkBkQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AkBQkDk")))  return HDGaussian_AkBQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AkBQkD")))   return HDGaussian_AkBQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AkBQDk")))   return HDGaussian_AkBQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AkBQD")))    return HDGaussian_AkBQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjBkQkDk"))) return HDGaussian_AjBkQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjBkQkD")))  return HDGaussian_AjBkQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjBkQDk")))  return HDGaussian_AjBkQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjBkQD")))   return HDGaussian_AjBkQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjBQkDk")))  return HDGaussian_AjBQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_AjBQkD")))   return HDGaussian_AjBQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_ABkQkDk")))  return HDGaussian_ABkQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_ABkQkD")))   return HDGaussian_ABkQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_ABkQDk")))   return HDGaussian_ABkQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_ABkQD")))    return HDGaussian_ABkQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_ABQkDk")))   return HDGaussian_ABQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_ABQkD")))    return HDGaussian_ABQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_ABQD")))     return HDGaussian_ABQD_;
  // try second naming sheme
  bool freeProp;
  Mixture res = stringToMixture(type, freeProp);
#ifdef STK_MIXTURE_DEBUG
  if (res == unknown_mixture_)
  { stk_cout << _T("In stringToMixture, mixture ") << type << _T(" not found.\n");}
#endif
  return res;
}
/* @ingroup Clustering
 *  convert a String to an Mixture.
 *  @param type the String we want to convert
 *  @return the Mixture represented by the String @c type. if the string
 *  does not match any known name, the @c unknown_mixture_ type is returned.
 **/
Mixture stringToMixture( String const& type, bool& freeProp)
{
  freeProp = false;
  if (toUpperString(type) == toUpperString(_T("Gamma_p_ajk_bjk"))) return Gamma_ajk_bjk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_p_ajk_bk")))  return Gamma_ajk_bk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_p_ajk_bj")))  return Gamma_ajk_bj_;
  if (toUpperString(type) == toUpperString(_T("Gamma_p_ajk_b")))   return Gamma_ajk_b_;
  if (toUpperString(type) == toUpperString(_T("Gamma_p_ak_bjk")))  return Gamma_ak_bjk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_p_ak_bk")))   return Gamma_ak_bk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_p_ak_bj")))   return Gamma_ak_bj_;
  if (toUpperString(type) == toUpperString(_T("Gamma_p_ak_b")))    return Gamma_ak_b_;
  if (toUpperString(type) == toUpperString(_T("Gamma_p_aj_bjk")))  return Gamma_aj_bjk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_p_aj_bk")))   return Gamma_aj_bk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_p_a_bjk")))   return Gamma_a_bjk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_p_a_bk")))    return Gamma_a_bk_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_p_sjk")))  return Gaussian_sjk_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_p_sk")))   return Gaussian_sk_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_p_sj")))   return Gaussian_sj_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_p_s")))    return Gaussian_s_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_p_sjsk"))) return Gaussian_sjsk_;
  if (toUpperString(type) == toUpperString(_T("Categorical_p_pjk"))) return Categorical_pjk_;
  if (toUpperString(type) == toUpperString(_T("Categorical_p_pk")))  return Categorical_pk_;
  if (toUpperString(type) == toUpperString(_T("Poisson_p_ljk")))  return Poisson_ljk_;
  if (toUpperString(type) == toUpperString(_T("Poisson_p_lk")))   return Poisson_lk_;
  if (toUpperString(type) == toUpperString(_T("Poisson_p_ljlk"))) return Poisson_ljlk_;
  if (toUpperString(type) == toUpperString(_T("Kmm_p_sk")))       return Kmm_sk_;
  if (toUpperString(type) == toUpperString(_T("Kmm_p_s")))        return Kmm_s_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjkBkQkDk"))) return HDGaussian_AjkBkQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjkBkQkD")))  return HDGaussian_AjkBkQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjkBkQDk")))  return HDGaussian_AjkBkQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjkBkQD")))   return HDGaussian_AjkBkQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjkBQkDk")))  return HDGaussian_AjkBQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjkBQkD")))   return HDGaussian_AjkBQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjkBQDk")))   return HDGaussian_AjkBQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjkBQD")))    return HDGaussian_AjkBQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AkBkQkDk")))  return HDGaussian_AkBkQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AkBkQkD")))   return HDGaussian_AkBkQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AkBkQDk")))   return HDGaussian_AkBkQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AkBkQD")))    return HDGaussian_AkBkQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AkBQkDk")))   return HDGaussian_AkBQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AkBQkD")))    return HDGaussian_AkBQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AkBQDk")))    return HDGaussian_AkBQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AkBQD")))     return HDGaussian_AkBQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjBkQkDk")))  return HDGaussian_AjBkQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjBkQkD")))   return HDGaussian_AjBkQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjBkQDk")))   return HDGaussian_AjBkQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjBkQD")))    return HDGaussian_AjBkQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjBQkDk")))   return HDGaussian_AjBQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_AjBQkD")))    return HDGaussian_AjBQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_ABkQkDk")))   return HDGaussian_ABkQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_ABkQkD")))    return HDGaussian_ABkQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_ABkQDk")))    return HDGaussian_ABkQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_ABkQD")))     return HDGaussian_ABkQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_ABQkDk")))    return HDGaussian_ABQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_ABQkD")))     return HDGaussian_ABQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_p_ABQD")))      return HDGaussian_ABQD_;
  freeProp = true;
  if (toUpperString(type) == toUpperString(_T("Gamma_pk_ajk_bjk"))) return Gamma_ajk_bjk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_pk_ajk_bk")))  return Gamma_ajk_bk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_pk_ajk_bj")))  return Gamma_ajk_bj_;
  if (toUpperString(type) == toUpperString(_T("Gamma_pk_ajk_b")))   return Gamma_ajk_b_;
  if (toUpperString(type) == toUpperString(_T("Gamma_pk_ak_bjk")))  return Gamma_ak_bjk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_pk_ak_bk")))   return Gamma_ak_bk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_pk_ak_bj")))   return Gamma_ak_bj_;
  if (toUpperString(type) == toUpperString(_T("Gamma_pk_ak_b")))    return Gamma_ak_b_;
  if (toUpperString(type) == toUpperString(_T("Gamma_pk_aj_bjk")))  return Gamma_aj_bjk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_pk_aj_bk")))   return Gamma_aj_bk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_pk_a_bjk")))   return Gamma_a_bjk_;
  if (toUpperString(type) == toUpperString(_T("Gamma_pk_a_bk")))    return Gamma_a_bk_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_pk_sjk")))  return Gaussian_sjk_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_pk_sk")))   return Gaussian_sk_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_pk_sj")))   return Gaussian_sj_;
  if (toUpperString(type) == toUpperString(_T("Gaussian_pk_s")))    return Gaussian_s_;
  if (toUpperString(type) == toUpperString(_T("Categorical_pk_pjk"))) return Categorical_pjk_;
  if (toUpperString(type) == toUpperString(_T("Categorical_pk_pk")))  return Categorical_pk_;
  if (toUpperString(type) == toUpperString(_T("Poisson_pk_ljk")))   return Poisson_ljk_;
  if (toUpperString(type) == toUpperString(_T("Poisson_pk_lk")))    return Poisson_lk_;
  if (toUpperString(type) == toUpperString(_T("Poisson_pk_ljlk")))  return Poisson_ljlk_;
  if (toUpperString(type) == toUpperString(_T("Kmm_pk_sk")))        return Kmm_sk_;
  if (toUpperString(type) == toUpperString(_T("Kmm_pk_s")))         return Kmm_s_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjkBkQkDk"))) return HDGaussian_AjkBkQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjkBkQkD")))  return HDGaussian_AjkBkQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjkBkQDk")))  return HDGaussian_AjkBkQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjkBkQD")))   return HDGaussian_AjkBkQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjkBQkDk")))  return HDGaussian_AjkBQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjkBQkD")))   return HDGaussian_AjkBQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjkBQDk")))   return HDGaussian_AjkBQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjkBQD")))    return HDGaussian_AjkBQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AkBkQkDk")))  return HDGaussian_AkBkQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AkBkQkD")))   return HDGaussian_AkBkQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AkBkQDk")))   return HDGaussian_AkBkQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AkBkQD")))    return HDGaussian_AkBkQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AkBQkDk")))   return HDGaussian_AkBQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AkBQkD")))    return HDGaussian_AkBQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AkBQDk")))    return HDGaussian_AkBQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AkBQD")))     return HDGaussian_AkBQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjBkQkDk")))  return HDGaussian_AjBkQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjBkQkD")))   return HDGaussian_AjBkQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjBkQDk")))   return HDGaussian_AjBkQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjBkQD")))    return HDGaussian_AjBkQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjBQkDk")))   return HDGaussian_AjBQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_AjBQkD")))    return HDGaussian_AjBQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_ABkQkDk")))   return HDGaussian_ABkQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_ABkQkD")))    return HDGaussian_ABkQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_ABkQDk")))    return HDGaussian_ABkQDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_ABkQD")))     return HDGaussian_ABkQD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_ABQkDk")))    return HDGaussian_ABQkDk_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_ABQkD")))     return HDGaussian_ABQkD_;
  if (toUpperString(type) == toUpperString(_T("HDGaussian_pk_ABQD")))      return HDGaussian_ABQD_;
#ifdef STK_MIXTURE_DEBUG
  stk_cout << _T("In stringToMixture, mixture ") << type << _T(" not found.\n");
#endif
  return unknown_mixture_;
}

/* @ingroup Clustering
 *  convert a Mixture to a String.
 *  @param type the type of Mixture we want to convert
 *  @return the string associated to this type.
 **/
String mixtureToString( Mixture const& type)
{
  if (type == Gamma_ajk_bjk_)     return String(_T("Gamma_ajk_bjk"));
  if (type == Gamma_ajk_bk_)      return String(_T("Gamma_ajk_bk"));
  if (type == Gamma_ajk_bj_)      return String(_T("Gamma_ajk_bj"));
  if (type == Gamma_ajk_b_)       return String(_T("Gamma_ajk_b"));
  if (type == Gamma_ak_bjk_)      return String(_T("Gamma_ak_bjk"));
  if (type == Gamma_ak_bk_)       return String(_T("Gamma_ak_bk"));
  if (type == Gamma_ak_bj_)       return String(_T("Gamma_ak_bj"));
  if (type == Gamma_ak_b_)        return String(_T("Gamma_ak_b"));
  if (type == Gamma_aj_bjk_)      return String(_T("Gamma_aj_bjk"));
  if (type == Gamma_aj_bk_)       return String(_T("Gamma_aj_bk"));
  if (type == Gamma_a_bjk_)       return String(_T("Gamma_a_bjk"));
  if (type == Gamma_a_bk_)        return String(_T("Gamma_a_bk"));
  if (type == Gaussian_sjk_)      return String(_T("Gaussian_sjk"));
  if (type == Gaussian_sk_)       return String(_T("Gaussian_sk"));
  if (type == Gaussian_sj_)       return String(_T("Gaussian_sj"));
  if (type == Gaussian_s_)        return String(_T("Gaussian_s"));
  if (type == Gaussian_sjsk_)     return String(_T("Gaussian_sjsk"));
  if (type == Categorical_pjk_)   return String(_T("Categorical_pjk"));
  if (type == Categorical_pk_)    return String(_T("Categorical_pk"));
  if (type == Poisson_ljk_)       return String(_T("Poisson_ljk"));
  if (type == Poisson_lk_)        return String(_T("Poisson_lk"));
  if (type == Poisson_ljlk_)      return String(_T("Poisson_ljlk"));
  if (type == Kmm_sk_)            return String(_T("Kmm_sk"));
  if (type == Kmm_s_)             return String(_T("Kmm_s"));
  if (type == HDGaussian_AjkBkQkDk_) return String(_T("HDGaussian_AjkBkQkDk"));
  if (type == HDGaussian_AjkBkQkD_)  return String(_T("HDGaussian_AjkBkQkD"));
  if (type == HDGaussian_AjkBkQDk_)  return String(_T("HDGaussian_AjkBkQDk"));
  if (type == HDGaussian_AjkBkQD_)   return String(_T("HDGaussian_AjkBkQD"));
  if (type == HDGaussian_AjkBQkDk_)  return String(_T("HDGaussian_AjkBQkDk"));
  if (type == HDGaussian_AjkBQkD_)   return String(_T("HDGaussian_AjkBQkD"));
  if (type == HDGaussian_AjkBQDk_)   return String(_T("HDGaussian_AjkBQDk"));
  if (type == HDGaussian_AjkBQD_)    return String(_T("HDGaussian_AjkBQD"));
  if (type == HDGaussian_AkBkQkDk_)  return String(_T("HDGaussian_AkBkQkDk"));
  if (type == HDGaussian_AkBkQkD_)   return String(_T("HDGaussian_AkBkQkD"));
  if (type == HDGaussian_AkBkQDk_)   return String(_T("HDGaussian_AkBkQDk"));
  if (type == HDGaussian_AkBkQD_)    return String(_T("HDGaussian_AkBkQD"));
  if (type == HDGaussian_AkBQkDk_)   return String(_T("HDGaussian_AkBQkDk"));
  if (type == HDGaussian_AkBQkD_)    return String(_T("HDGaussian_AkBQkD"));
  if (type == HDGaussian_AkBQDk_)    return String(_T("HDGaussian_AkBQDk"));
  if (type == HDGaussian_AkBQD_)     return String(_T("HDGaussian_AkBQD"));
  if (type == HDGaussian_AjBkQkDk_)  return String(_T("HDGaussian_AjBkQkDk"));
  if (type == HDGaussian_AjBkQkD_)   return String(_T("HDGaussian_AjBkQkD"));
  if (type == HDGaussian_AjBkQDk_)   return String(_T("HDGaussian_AjBkQDk"));
  if (type == HDGaussian_AjBkQD_)    return String(_T("HDGaussian_AjBkQD"));
  if (type == HDGaussian_AjBQkDk_)   return String(_T("HDGaussian_AjBQkDk"));
  if (type == HDGaussian_AjBQkD_)    return String(_T("HDGaussian_AjBQkD"));
  if (type == HDGaussian_ABkQkDk_)   return String(_T("HDGaussian_ABkQkDk"));
  if (type == HDGaussian_ABkQkD_)    return String(_T("HDGaussian_ABkQkD"));
  if (type == HDGaussian_ABkQDk_)    return String(_T("HDGaussian_ABkQDk"));
  if (type == HDGaussian_ABkQD_)     return String(_T("HDGaussian_ABkQD"));
  if (type == HDGaussian_ABQkDk_)    return String(_T("HDGaussian_ABQkDk"));
  if (type == HDGaussian_ABQkD_)     return String(_T("HDGaussian_ABQkD"));
  if (type == HDGaussian_ABQD_)      return String(_T("HDGaussian_ABQD"));
  return String(_T("unknown"));
}

/* @ingroup Clustering
 *  convert a Mixture to a string specifying if the model is with free
 *  proportions.
 *  @sa stringToMixture
 *  @param type the Mixture we want to convert
 *  @param freeProp @c true if the model have free proportions, @c false otherwise.
 *  @return the string represented by the Mixture @c type.
 **/
String mixtureToString(Mixture type, bool freeProp)
{
  if (freeProp == false)
  {
    if (type == Gamma_ajk_bjk_)     return String(_T("Gamma_p_ajk_bjk"));
    if (type == Gamma_ajk_bk_)      return String(_T("Gamma_p_ajk_bk"));
    if (type == Gamma_ajk_bj_)      return String(_T("Gamma_p_ajk_bj"));
    if (type == Gamma_ajk_b_)       return String(_T("Gamma_p_ajk_b"));
    if (type == Gamma_ak_bjk_)      return String(_T("Gamma_p_ak_bjk"));
    if (type == Gamma_ak_bk_)       return String(_T("Gamma_p_ak_bk"));
    if (type == Gamma_ak_bj_)       return String(_T("Gamma_p_ak_bj"));
    if (type == Gamma_ak_b_)        return String(_T("Gamma_p_ak_b"));
    if (type == Gamma_aj_bjk_)      return String(_T("Gamma_p_aj_bjk"));
    if (type == Gamma_aj_bk_)       return String(_T("Gamma_p_aj_bk"));
    if (type == Gamma_a_bjk_)       return String(_T("Gamma_p_a_bjk"));
    if (type == Gamma_a_bk_)        return String(_T("Gamma_p_a_bk"));
    if (type == Gaussian_sjk_)      return String(_T("Gaussian_p_sjk"));
    if (type == Gaussian_sk_)       return String(_T("Gaussian_p_sk"));
    if (type == Gaussian_sj_)       return String(_T("Gaussian_p_sj"));
    if (type == Gaussian_s_)        return String(_T("Gaussian_p_s"));
    if (type == Gaussian_sjsk_)     return String(_T("Gaussian_p_sjsk"));
    if (type == Categorical_pjk_)   return String(_T("Categorical_p_pjk"));
    if (type == Categorical_pk_)    return String(_T("Categorical_p_pk"));
    if (type == Poisson_ljk_)       return String(_T("Poisson_p_ljk"));
    if (type == Poisson_lk_)        return String(_T("Poisson_p_lk"));
    if (type == Poisson_ljlk_)      return String(_T("Poisson_p_ljlk"));
    if (type == Poisson_ljlk_)      return String(_T("Poisson_p_ljlk"));
    if (type == Kmm_sk_)            return String(_T("Kmm_p_sk"));
    if (type == Kmm_s_)             return String(_T("Kmm_p_s"));
    if (type == HDGaussian_AjkBkQkDk_) return String(_T("HDGaussian_p_AjkBkQkDk"));
    if (type == HDGaussian_AjkBkQkD_)  return String(_T("HDGaussian_p_AjkBkQkD"));
    if (type == HDGaussian_AjkBkQDk_)  return String(_T("HDGaussian_p_AjkBkQDk"));
    if (type == HDGaussian_AjkBkQD_)   return String(_T("HDGaussian_p_AjkBkQD"));
    if (type == HDGaussian_AjkBQkDk_)  return String(_T("HDGaussian_p_AjkBQkDk"));
    if (type == HDGaussian_AjkBQkD_)   return String(_T("HDGaussian_p_AjkBQkD"));
    if (type == HDGaussian_AjkBQDk_)   return String(_T("HDGaussian_p_AjkBQDk"));
    if (type == HDGaussian_AjkBQD_)    return String(_T("HDGaussian_p_AjkBQD"));
    if (type == HDGaussian_AkBkQkDk_)  return String(_T("HDGaussian_p_AkBkQkDk"));
    if (type == HDGaussian_AkBkQkD_)   return String(_T("HDGaussian_p_AkBkQkD"));
    if (type == HDGaussian_AkBkQDk_)   return String(_T("HDGaussian_p_AkBkQDk"));
    if (type == HDGaussian_AkBkQD_)    return String(_T("HDGaussian_p_AkBkQD"));
    if (type == HDGaussian_AkBQkDk_)   return String(_T("HDGaussian_p_AkBQkDk"));
    if (type == HDGaussian_AkBQkD_)    return String(_T("HDGaussian_p_AkBQkD"));
    if (type == HDGaussian_AkBQDk_)    return String(_T("HDGaussian_p_AkBQDk"));
    if (type == HDGaussian_AkBQD_)     return String(_T("HDGaussian_p_AkBQD"));
    if (type == HDGaussian_AjBkQkDk_)  return String(_T("HDGaussian_p_AjBkQkDk"));
    if (type == HDGaussian_AjBkQkD_)   return String(_T("HDGaussian_p_AjBkQkD"));
    if (type == HDGaussian_AjBkQDk_)   return String(_T("HDGaussian_p_AjBkQDk"));
    if (type == HDGaussian_AjBkQD_)    return String(_T("HDGaussian_p_AjBkQD"));
    if (type == HDGaussian_AjBQkDk_)   return String(_T("HDGaussian_p_AjBQkDk"));
    if (type == HDGaussian_AjBQkD_)    return String(_T("HDGaussian_p_AjBQkD"));
    if (type == HDGaussian_ABkQkDk_)   return String(_T("HDGaussian_p_ABkQkDk"));
    if (type == HDGaussian_ABkQkD_)    return String(_T("HDGaussian_p_ABkQkD"));
    if (type == HDGaussian_ABkQDk_)    return String(_T("HDGaussian_p_ABkQDk"));
    if (type == HDGaussian_ABkQD_)     return String(_T("HDGaussian_p_ABkQD"));
    if (type == HDGaussian_ABQkDk_)    return String(_T("HDGaussian_p_ABQkDk"));
    if (type == HDGaussian_ABQkD_)     return String(_T("HDGaussian_p_ABQkD"));
    if (type == HDGaussian_ABQD_)     return String(_T("HDGaussian_p_ABQD"));
  }
  else
  {
    if (type == Gamma_ajk_bjk_)     return String(_T("Gamma_pk_ajk_bjk"));
    if (type == Gamma_ajk_bk_)      return String(_T("Gamma_pk_ajk_bk"));
    if (type == Gamma_ajk_bj_)      return String(_T("Gamma_pk_ajk_bj"));
    if (type == Gamma_ajk_b_)       return String(_T("Gamma_pk_ajk_b"));
    if (type == Gamma_ak_bjk_)      return String(_T("Gamma_pk_ak_bjk"));
    if (type == Gamma_ak_bk_)       return String(_T("Gamma_pk_ak_bk"));
    if (type == Gamma_ak_bj_)       return String(_T("Gamma_pk_ak_bj"));
    if (type == Gamma_ak_b_)        return String(_T("Gamma_pk_ak_b"));
    if (type == Gamma_aj_bjk_)      return String(_T("Gamma_pk_aj_bjk"));
    if (type == Gamma_aj_bk_)       return String(_T("Gamma_pk_aj_bk"));
    if (type == Gamma_a_bjk_)       return String(_T("Gamma_pk_a_bjk"));
    if (type == Gamma_a_bk_)        return String(_T("Gamma_pk_a_bk"));
    if (type == Gaussian_sjk_)      return String(_T("Gaussian_pk_sjk"));
    if (type == Gaussian_sk_)       return String(_T("Gaussian_pk_sk"));
    if (type == Gaussian_sj_)       return String(_T("Gaussian_pk_sj"));
    if (type == Gaussian_s_)        return String(_T("Gaussian_pk_s"));
    if (type == Gaussian_sjsk_)     return String(_T("Gaussian_pk_sjsk"));
    if (type == Categorical_pjk_)   return String(_T("Categorical_pk_pjk"));
    if (type == Categorical_pk_)    return String(_T("Categorical_pk_pk"));
    if (type == Poisson_ljk_)       return String(_T("Poisson_pk_ljk"));
    if (type == Poisson_lk_)        return String(_T("Poisson_pk_lk"));
    if (type == Poisson_ljlk_)      return String(_T("Poisson_pk_ljlk"));
    if (type == Kmm_sk_)            return String(_T("Kmm_pk_sk"));
    if (type == Kmm_s_)             return String(_T("Kmm_pk_s"));
    if (type == HDGaussian_AjkBkQkDk_) return String(_T("HDGaussian_pk_AjkBkQkDk"));
    if (type == HDGaussian_AjkBkQkD_)  return String(_T("HDGaussian_pk_AjkBkQkD"));
    if (type == HDGaussian_AjkBkQDk_)  return String(_T("HDGaussian_pk_AjkBkQDk"));
    if (type == HDGaussian_AjkBkQD_)   return String(_T("HDGaussian_pk_AjkBkQD"));
    if (type == HDGaussian_AjkBQkDk_)  return String(_T("HDGaussian_pk_AjkBQkDk"));
    if (type == HDGaussian_AjkBQkD_)   return String(_T("HDGaussian_pk_AjkBQkD"));
    if (type == HDGaussian_AjkBQDk_)   return String(_T("HDGaussian_pk_AjkBQDk"));
    if (type == HDGaussian_AjkBQD_)    return String(_T("HDGaussian_pk_AjkBQD"));
    if (type == HDGaussian_AkBkQkDk_)  return String(_T("HDGaussian_pk_AkBkQkDk"));
    if (type == HDGaussian_AkBkQkD_)   return String(_T("HDGaussian_pk_AkBkQkD"));
    if (type == HDGaussian_AkBkQDk_)   return String(_T("HDGaussian_pk_AkBkQDk"));
    if (type == HDGaussian_AkBkQD_)    return String(_T("HDGaussian_pk_AkBkQD"));
    if (type == HDGaussian_AkBQkDk_)   return String(_T("HDGaussian_pk_AkBQkDk"));
    if (type == HDGaussian_AkBQkD_)    return String(_T("HDGaussian_pk_AkBQkD"));
    if (type == HDGaussian_AkBQDk_)    return String(_T("HDGaussian_pk_AkBQDk"));
    if (type == HDGaussian_AkBQD_)     return String(_T("HDGaussian_pk_AkBQD"));
    if (type == HDGaussian_AjBkQkDk_)  return String(_T("HDGaussian_pk_AjBkQkDk"));
    if (type == HDGaussian_AjBkQkD_)   return String(_T("HDGaussian_pk_AjBkQkD"));
    if (type == HDGaussian_AjBkQDk_)   return String(_T("HDGaussian_pk_AjBkQDk"));
    if (type == HDGaussian_AjBkQD_)    return String(_T("HDGaussian_pk_AjBkQD"));
    if (type == HDGaussian_AjBQkDk_)   return String(_T("HDGaussian_pk_AjBQkDk"));
    if (type == HDGaussian_AjBQkD_)    return String(_T("HDGaussian_pk_AjBQkD"));
    if (type == HDGaussian_ABkQkDk_)   return String(_T("HDGaussian_pk_ABkQkDk"));
    if (type == HDGaussian_ABkQkD_)    return String(_T("HDGaussian_pk_ABkQkD"));
    if (type == HDGaussian_ABkQDk_)    return String(_T("HDGaussian_pk_ABkQDk"));
    if (type == HDGaussian_ABkQD_)     return String(_T("HDGaussian_pk_ABkQD"));
    if (type == HDGaussian_ABQkDk_)    return String(_T("HDGaussian_pk_ABQkDk"));
    if (type == HDGaussian_ABQkD_)     return String(_T("HDGaussian_pk_ABQkD"));
    if (type == HDGaussian_ABQD_)      return String(_T("HDGaussian_pk_ABQD"));
  }
  return String(_T("unknown"));
}

IMixtureCriterion* createCriterion( Clust::criterionType type)
{
  IMixtureCriterion* p_criter = 0;
  switch (type)
  {
  case aic_:
    p_criter = new AICMixtureCriterion();
    break;
  case bic_:
    p_criter = new BICMixtureCriterion();
    break;
  case icl_:
    p_criter = new ICLMixtureCriterion();
    break;
  case ml_:
    p_criter = new MLMixtureCriterion();
    break;
  }
  return p_criter;
}

/* @return a pointer on the class computing the criterion */
STK::IMixtureCriterion* createCriterion( String const& criterion)
{
  criterionType type = STK::Clust::stringToCriterion(criterion);
  return createCriterion(type);
}

/* utility function for creating an estimation algorithm
 *  @param algo the algorithm to create
 *  @param nbIterMax the maximal number of iteration of the algorithm
 *  @param epsilon the tolerance of the algorithm
 **/
IMixtureAlgo* createAlgo(Clust::algoType algo, int nbIterMax, Real epsilon)
{
  IMixtureAlgo* p_algo = 0;
  switch (algo)
  {
  case emAlgo_:
    p_algo = new EMAlgo();
    break;
  case cemAlgo_:
    p_algo = new CEMAlgo();
    break;
  case semAlgo_:
    p_algo = new SEMAlgo();
    break;
  case semiSemAlgo_:
    p_algo = new SemiSEMAlgo();
    break;
  default:
    break;
  }
  if (p_algo)
  {
    p_algo->setNbIterMax(nbIterMax);
    p_algo->setEpsilon(epsilon);
  }
  return p_algo;
}

/* utility function for creating an estimation algorithm
 *  @param algo the algorithm to create
 *  @param nbIterMax the maximal number of iteration of the algorithm
 *  @param epsilon the tolerance of the algorithm
 **/
IMixtureAlgoLearn* createLearnAlgo(Clust::algoLearnType algo, int nbIterMax, Real epsilon)
{
  IMixtureAlgoLearn* p_algo = 0;
  switch (algo)
  {
  case imputeAlgo_:
    p_algo = new ImputeAlgo();
    break;
  case simulAlgo_:
    p_algo = new SimulAlgo();
    break;
  default:
    break;
  }
  if (p_algo)
  {
    p_algo->setNbIterMax(nbIterMax);
    p_algo->setEpsilon(epsilon);
  }
  return p_algo;
}

/* utility function for creating an estimation algorithm
 *  @param algo the algorithm to create
 *  @param nbIterMax the maximal number of iteration of the algorithm
 *  @param epsilon the tolerance of the algorithm
 **/
IMixtureAlgoPredict* createPredictAlgo(Clust::algoPredictType algo, int nbIterBurn, int nbIterLong, Real epsilon)
{
  IMixtureAlgoPredict* p_algo = 0;
  switch (algo)
  {
    case emPredictAlgo_:
      p_algo = new EMPredict();
      break;
    case semiSEMPredictAlgo_:
      p_algo = new SemiSEMPredict();
      break;
    default:
      break;
  }
  if (p_algo)
  {
    p_algo->setNbIterBurn(nbIterBurn);
    p_algo->setNbIterLong(nbIterLong);
    p_algo->setEpsilon(epsilon);
  }
  return p_algo;
}

/* utility function for creating a model initializer
 *  @param init the kind of initializer to create
 *  @param algo the kind of algorithm to add to the initializer
 *  @param nbIterMax the maximal number of iteration of the algorithm
 *  @param epsilon the tolerance of the algorithm
 **/
IMixtureInit* createInit(Clust::initType init,  int nbInits, Clust::algoType algo, int nbIterMax, Real epsilon)
{
  IMixtureInit* p_init = 0;
  switch (init)
  {
    case Clust::noInit_:
      p_init = 0;
      break;
    case Clust::randomInit_:
      p_init = new RandomInit();
      break;
    case Clust::randomParamInit_:
      p_init = new RandomInit();
      break;
    case Clust::randomClassInit_:
      p_init = new ClassInit();
      break;
    case Clust::randomFuzzyInit_:
      p_init = new FuzzyInit();
      break;
    default:
      break;
  }
  if (p_init)
  {
    p_init->setNbTry(nbInits);
    p_init->setInitAlgo(Clust::createAlgo(algo, nbIterMax, epsilon));
  }
  return p_init;
}

/* @ingroup Clustering
 *  Utility function for creating a SimpleStrategy.
 *  @param nbTry the number of tries.
 *  @param p_init the initializer to use.
 *  @param algo the algorithm to use in the long run.
 *  @return an instance of the SimpleStrategy
 **/
IMixtureStrategy* createSimpleStrategy( IMixtureComposer*& p_composer
                                      , int nbTry
                                      , IMixtureInit* const& p_init
                                      , IMixtureAlgo* const& algo)
{
  SimpleStrategyParam* p_strategyParam = new SimpleStrategyParam();
  p_strategyParam->p_algo_ = algo;
  SimpleStrategy* p_strategy = new SimpleStrategy(p_composer);
  p_strategy->setNbTry(nbTry);
  p_strategy->setMixtureInit(p_init);
  p_strategy->setParam(p_strategyParam);
  return p_strategy;
}


/* @ingroup Clustering
 *  Utility function for creating a FullStrategy.
 *  @param nbTry the number of tries.
 *  @param p_init the initializer to use.
 *  @param nbShortRun the number of shortRun.
 *  @param shortRunAlgo the algorithm to use in the short run.
 *  @param longRunAlgo the algorithm to use in the long run.
 *  @return an instance of the SimpleStrategy
 **/
IMixtureStrategy* createFullStrategy( IMixtureComposer*& p_composer
                                    , int nbTry, int nbInitRun
                                    , IMixtureInit* const& p_init
                                    , int nbShortRun, IMixtureAlgo* const& shortRunAlgo
                                    , IMixtureAlgo* const& longRunAlgo)
{
  FullStrategyParam* p_strategyParam = new FullStrategyParam();
  p_strategyParam->nbInitRun_  = nbInitRun;
  p_strategyParam->nbShortRun_  = nbShortRun;
  p_strategyParam->p_shortAlgo_ = shortRunAlgo;
  p_strategyParam->p_longAlgo_  = longRunAlgo;
  FullStrategy* p_strategy = new FullStrategy(p_composer);
  p_strategy->setNbTry(nbTry);
  p_strategy->setMixtureInit(p_init);
  p_strategy->setParam(p_strategyParam);
  return p_strategy;
}

} // namespace Clust

} // namespace STK


