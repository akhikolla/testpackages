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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::STatistiK::Law
 * created on: 2 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Law_Util.h
 *  @brief In this file we define the utilities constant and method for the Law
 *  namespace.
 **/


#ifndef STK_LAW_UTIL_H
#define STK_LAW_UTIL_H

#include "STK_RandBase.h"

namespace STK
{

namespace Law
{

/** default random number generator. */
static RandBase generator;

/** @ingroup Law
 * list of the univariate distribution laws
 **/
enum UnivariateDistribution
{
  bernoulli_,
  beta_,
  binomial_,
  categorical_,
  cauchy_,
  chisquared_,
  exponential_,
  fisher_snedecor_,
  gamma_,
  geometric_,
  hypergeometric_,
  logistic_,
  lognormal_,
  negative_binomial_,
  normal_,
  poisson_,
  student_,
  uniform_,
  uniform_discrete_,
  weibull_,
  unknown_univ_distribution_
};


/** @ingroup Law
 *  Convert a String to a UnivariateDistribution. The recognized strings are
 * <table >
 * <tr> <th> Distribution     </th> </tr>
 * <tr> <td> "beta"           </td></tr>
 * <tr> <td> "binomial"       </td></tr>
 * <tr> <td> "categorical"    </td></tr>
 * <tr> <td> "cauchy"         </td></tr>
 * <tr> <td> "chi squared"    </td></tr>
 * <tr> <td> "exponential"    </td></tr>
 * <tr> <td> "fisher snedecor" </td></tr>
 * <tr> <td> "gamma"          </td></tr>
 * <tr> <td> "geometric"      </td></tr>
 * <tr> <td> "hypergeometric" </td></tr>
 * <tr> <td> "logistic"       </td></tr>
 * <tr> <td> "lognormal"      </td></tr>
 * <tr> <td> "negative binomial" </td></tr>
 * <tr> <td> "normal"         </td></tr>
 * <tr> <td> "poisson"        </td></tr>
 * <tr> <td> "student"        </td></tr>
 * <tr> <td> "uniform"        </td></tr>
 * <tr> <td> "uniform discrete" </td></tr>
 * <tr> <td> "weibull"          </td></tr>
 * </table>
 *  @param dist the String we want to convert
 *  @return the Mixture represented by the String @c type. if the string
 *  does not match any known name, the @c unknown_mixture_ type is returned.
 **/
inline UnivariateDistribution stringToUnivariateDistribution( std::string const& dist)
{
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("beta"))) return beta_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("binomial"))) return binomial_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("categorical"))) return categorical_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("cauchy"))) return cauchy_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("chisquared"))) return chisquared_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("exponential"))) return exponential_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("fishersnedecor"))) return fisher_snedecor_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("gamma"))) return gamma_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("geometric"))) return geometric_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("hypergeometric"))) return hypergeometric_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("logistic"))) return logistic_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("lognormal"))) return lognormal_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("negativebinomial"))) return negative_binomial_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("normal"))) return normal_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("poisson"))) return poisson_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("student"))) return student_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("uniform"))) return uniform_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("uniformdiscrete"))) return uniform_discrete_;
  if (toUpperString(removeWhiteSpaces(dist)) == toUpperString(_T("weibull"))) return weibull_;
#ifdef STK_STATISTIK_DEBUG
  stk_cout << _T("In stringToUnivariateDistribution, dist ") << dist << _T(" not found.\n");
#endif
  return unknown_univ_distribution_;
}

/** @ingroup Law
 *  convert a UnivariateDistribution to a String.
 *  @param type the type of UnivariateDistribution we want to convert
 *  @return the string associated to this type.
 *  @sa stringToUnivariateDistribution
 **/
inline std::string univariateDistributionToString( UnivariateDistribution const& type)
{
  if (type == beta_) return String(_T("beta"));
  if (type == binomial_) return String(_T("binomial"));
  if (type == categorical_) return String(_T("categorical"));
  if (type == chisquared_) return String(_T("chi squared"));
  if (type == cauchy_) return String(_T("cauchy"));
  if (type == exponential_) return String(_T("exponential"));
  if (type == fisher_snedecor_) return String(_T("fisher snedecor"));
  if (type == gamma_) return String(_T("gamma"));
  if (type == geometric_) return String(_T("geometric"));
  if (type == hypergeometric_) return String(_T("hypergeometric"));
  if (type == logistic_) return String(_T("logistic"));
  if (type == lognormal_) return String(_T("lognormal"));
  if (type == negative_binomial_) return String(_T("negative binomial"));
  if (type == normal_) return String(_T("normal"));
  if (type == poisson_) return String(_T("poisson"));
  if (type == student_) return String(_T("student"));
  if (type == uniform_) return String(_T("uniform"));
  if (type == uniform_discrete_) return String(_T("uniform discrete"));
  if (type == weibull_) return String(_T("weibull"));
  return String(_T("unknown"));
}



}  // namespace Law

}  // namespace STK

#endif /* STK_LAW_UTIL_H */
