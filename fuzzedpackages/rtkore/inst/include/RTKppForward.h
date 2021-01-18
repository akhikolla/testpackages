/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Inria

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
 * Project:  stkpp
 * created on: 22 août 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file RTKppForward.h
 *  @brief In this file we include the header files needed for the integration of stkpp to Rcpp.
 **/


#ifndef RTKPPFORWARD_H
#define RTKPPFORWARD_H

#include <RcppCommon.h>
#include <Rmath.h>
#include "STKpp.h"

/* forward declarations */
namespace STK
{
template <typename Type> class RVector;
template <typename Type> class RMatrix;

} // namespace STK


/* Rcpp integration */
namespace Rcpp
{
  /* This support for wrap, will not work as base template
   * template <typename T>
   * inline SEXP wrap(const T& object);
   * of Rcpp will be preferred by the compiler...
   *
   * If you want to wrap an expression, you have to use
   * STK::wrap()
   **/
  template<typename Derived>
  SEXP wrap( STK::ExprBase<Derived> const& matrix);

  template<typename Type>
  SEXP wrap( STK::Array2D<Type> const& matrix);
  template<typename Type>
  SEXP wrap( STK::Array2DVector<Type> const& vec);
  template<typename Type>
  SEXP wrap( STK::Array2DPoint<Type> const& vec);
  template <typename Type, int SizeRows_, int SizeCols_, bool Orient_>
  SEXP wrap( STK::CArray<Type, SizeRows_, SizeCols_, Orient_> const& matrix);
  template <typename Type, int SizeRows_, bool Orient_>
  SEXP wrap( STK::CArrayVector<Type, SizeRows_, Orient_ > const& vec);
  template <typename Type, int SizeCols_, bool Orient_>
  SEXP wrap( STK::CArrayPoint<Type, SizeCols_, Orient_> const& vec);

namespace traits
{
  /* support for as */
  template<typename Type> class Exporter< STK::RVector<Type> >;
  template<typename Type> class Exporter< STK::RMatrix<Type> >;
} // namespace traits

} // namespace Rcpp

#endif /* RTKPPFORWARD_H */
