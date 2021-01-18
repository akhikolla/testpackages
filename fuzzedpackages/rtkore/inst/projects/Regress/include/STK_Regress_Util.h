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
 * Project:  stkpp::Regress
 * created on: 23 juin 2011
 * Purpose:  num and other utilities for the Regress project..
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Regress_Util.h
 *  @brief In this file we declare the utilities enumerations and methods for
 *  the Regress project.
 **/


#ifndef STK_REGRESS_UTIL_H
#define STK_REGRESS_UTIL_H

#include <STKernel/include/STK_String.h>

namespace STK
{

namespace Regress
{

/** @ingroup Regress
 *  functional basis coefficients to use  for fonctional and non-linear
 *  regression
 **/
enum TypeBasisFunction
{
  /** unknown coefficients */
  unknown_basis_ = -1
  /** BSpline coefficients */
  , bspline_
  /** Sines coefficients */
  , sines_
  /** cosines coefficients */
  , cosines_
  /** trigonometric coefficients (sines and cosines) */
  , trigonometric_
  /**Chebyshev coefficients */
  , chebyshev
};

/** @ingroup Regress
 * Regression method to use. */
enum TypeRegression
{
  /** unknown regression */
  unknown_regression_ = -1
  /** linear regression */
  , linear_
  /** additive BSpline regression */
  , additiveBSpline_
  /** adaptive BSpline regression */
  , adaptiveBSpline_
};

/** @ingroup Regress
 *  Convert a String to a TypeRegression.
 *  @param type the String we want to convert
 *  @return the TypeRegression represented by the String @c type. if the string
 *  does not match any known name, the @c unknown_regression_ type is returned.
 **/
TypeRegression stringToTypeRegression( String const& type);

/** @ingroup Regress
 *  Convert a TypeRegression to a String.
 *  @param type the type of regression we want to convert
 *  @return the string associated to this type.
 **/
String typeRegressionToString( TypeRegression const& type);

/** @ingroup Regress
 *  Method to use for positioning the knots for BSpline basis
 **/
enum KnotsPosition
{
  uniformKnotsPositions_  ///< uniform knots
 , periodicKnotsPositions_ ///< periodic knots
 , densityKnotsPositions_  ///< knots using density of the data
 , unknown_Knots_Position_ = -1  ///< unknown positions
};

/** convert a String to a KnotsPosition.
 *  @param type the type of KnotsPosition in a string
 *  @return the KnotsPosition represented by the String @c type. If the string
 *  does not match any known name, the @c unknown_regression_ type is returned.
 **/
KnotsPosition stringToKnotsPosition( String const& type);

/** @ingroup Regress
 *  convert a KnotsPosition to a String.
 *  @param type the KnotsPosition we want to convert to a string
 *  @return the string associated to this KnotsPosition
 **/
String knotsPositionToString( KnotsPosition const& type);



//
// implementation
/* convert a String to a TypeReduction.
 *  @param type the type of reduction we want to define
 *  @return the TypeReduction represented by the String @c type. if the string
 *  does not match any known name, the @c unknown_regression_ type is returned.
 **/
inline KnotsPosition stringToKnotsPosition( String const& type)
{
  if (toUpperString(type) == toUpperString(_T("uniform")))  return uniformKnotsPositions_;
  if (toUpperString(type) == toUpperString(_T("periodic"))) return periodicKnotsPositions_;
  if (toUpperString(type) == toUpperString(_T("density"))) return densityKnotsPositions_;
  return unknown_Knots_Position_;
}

/* convert a TypeReduction to a String.
 *  @param type the type of reduction we want to convert
 *  @return the string associated to this type.
 **/
inline String knotsPositionToString( KnotsPosition const& type)
{
  if (type == uniformKnotsPositions_)  return String(_T("uniform"));
  if (type == periodicKnotsPositions_) return String(_T("periodic"));
  if (type == densityKnotsPositions_) return String(_T("density"));
  return String(_T("unknown"));
}

/* convert a String to a TypeRegression.
 *  @param type the type of regression we want to define
 *  @return the TypeRegression represented by the String @c type. if the string
 *  does not match any known name, the @c unknown_regression_ type is returned.
 **/
inline TypeRegression stringToTypeRegression( String const& type)
{
  if (toUpperString(type) == toUpperString(_T("unknown"))) return unknown_regression_;
  if (toUpperString(type) == toUpperString(_T("linear")))  return linear_;
  if (toUpperString(type) == toUpperString(_T("additiveBSpline"))) return additiveBSpline_;
  if (toUpperString(type) == toUpperString(_T("adaptiveBSpline"))) return adaptiveBSpline_;
  return unknown_regression_;
}

/* convert a TypeRegression to a String.
 *  @param type the type of regression we want to convert
 *  @return the string associated to this type.
 **/
inline String typeRegressionToString( TypeRegression const& type)
{
  if (type == linear_)  return String(_T("linear"));
  if (type == additiveBSpline_) return String(_T("additiveBSpline"));
  if (type == adaptiveBSpline_) return String(_T("adaptiveBSpline"));
  return String(_T("unknown"));
}

} // namespace Regress

} //namespace STK

#endif /* STK_REGRESS_UTIL_H */
