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

/** @file STK_Basis_Util.h
 *  @brief In this file we declare the utilities enumerations and methods for
 *  the Basis functions
 **/


#ifndef STK_BASIS_UTIL_H
#define STK_BASIS_UTIL_H

#include <STKernel/include/STK_String.h>

namespace STK
{

namespace Basis
{

/** @ingroup Regress
 *  function basis function to use for functional and non-linear regression
 **/
enum TypeBasisFunction
{
  /** BSpline basis */
  bspline_
  /** Sines basis */
  , sines_
  /** cosines basis */
  , cosines_
  /** trigonometric basis (sines and cosines) */
  , trigonometric_
  /**Chebyshev basis */
  , chebyshev_
  /** unknown basis */
  , unknown_basis_ = -1
};


/** @ingroup Regress
 *  Method to use for positioning the knots in BSpline basis */
enum KnotsPosition
{
   uniformKnotsPositions_         ///< uniform knots
 , periodicKnotsPositions_        ///< periodic knots
 , densityKnotsPositions_         ///< knots using density of the data
 , unknown_Knots_Position_ = -1   ///< unknown positions
};

/** @ingroup Regress
 *  convert a String to a KnotsPosition.
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

/** @ingroup Regress
 *  convert a String to a TypeBasisFunction.
 *  @param type the type of TypeBasisFunction in a string
 *  @return the TypeBasisFunction represented by the String @c type. If the string
 *  does not match any known name, the @c unknown_basis_ type is returned.
 **/
TypeBasisFunction stringToTypeBasisFunction( String const& type);

/** @ingroup Regress
 *  convert a TypeBasisFunction to a String.
 *  @param type the TypeBasisFunction we want to convert to a string
 *  @return the string associated to this TypeBasisFunction
 **/
String typeBasisFunctionToString( TypeBasisFunction const& type);


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
  if (toUpperString(type) == toUpperString(_T("density")))  return densityKnotsPositions_;
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

/* @ingroup Regress
 *  convert a String to a TypeBasisFunction.
 *  @param type the type of TypeBasisFunction in a string
 *  @return the TypeBasisFunction represented by the String @c type. If the string
 *  does not match any known name, the @c unknown_basis_ type is returned.
 **/
inline TypeBasisFunction stringToTypeBasisFunction( String const& type)
{
  if (toUpperString(type) == toUpperString(_T("bspline")))       return bspline_;
  if (toUpperString(type) == toUpperString(_T("sines")))         return sines_;
  if (toUpperString(type) == toUpperString(_T("cosines")))       return cosines_;
  if (toUpperString(type) == toUpperString(_T("trigonometric"))) return trigonometric_;
  if (toUpperString(type) == toUpperString(_T("chebyshev")))     return chebyshev_;
  return unknown_basis_;
}


/* @ingroup Regress
 *  convert a TypeBasisFunction to a String.
 *  @param type the TypeBasisFunction we want to convert to a string
 *  @return the string associated to this TypeBasisFunction
 **/
inline String typeBasisFunctionToString( TypeBasisFunction const& type)
{
  if (type == bspline_)       return String(_T("bspline"));
  if (type == sines_)         return String(_T("sines"));
  if (type == cosines_)       return String(_T("cosines"));
  if (type == trigonometric_) return String(_T("trigonometric"));
  if (type == chebyshev_)     return String(_T("chebyshev"));
  return String(_T("unknown"));
}


} // namespace Basis

} //namespace STK

#endif /* STK_BASIS_UTIL_H */
