/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

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

/** @file Rcpp_wrap.h
 *  @brief In this file we partially specialize the Rcpp::wrap converter.
 **/

#ifndef RCPP_WRAP_H
#define RCPP_WRAP_H

/* Rcpp integration */
namespace Rcpp
{
  /* support for wrap */
  template<typename Type>
  SEXP wrap( STK::Array2D<Type> const& matrix)
  {
    enum
    {
      Rtype_ = STK::hidden::RcppTraits<Type>::Rtype_
    };

    Matrix<Rtype_> res(matrix.sizeRows(), matrix.sizeCols());
    for(int j=matrix.beginCols(), jRes=0; j< matrix.endCols(); j++, jRes++)
    {
      for(int i=matrix.beginRows(), iRes=0; i< matrix.endRows(); i++, iRes++)
      { res(iRes, jRes) = matrix.elt(i,j);}
    }
    return Rcpp::wrap(res);
  }
  /* support for wrap */
  template<typename Type>
  SEXP wrap( STK::Array2DVector<Type> const& vec)
  {
    enum
    {
      Rtype_ = STK::hidden::RcppTraits<Type>::Rtype_
    };

    Vector<Rtype_> res(vec.size());
    for(int i=vec.begin(), iRes=0; i< vec.end(); i++, iRes++)
    { res(iRes) = vec.elt(i);}
    return Rcpp::wrap(res);
  }
  /* support for wrap */
  template<typename Type>
  SEXP wrap( STK::Array2DPoint<Type> const& vec)
  {
    enum
    {
      Rtype_ = STK::hidden::RcppTraits<Type>::Rtype_
    };
    Vector<Rtype_> res(vec.size());
    for(int i=vec.begin(), iRes=0; i< vec.end(); i++, iRes++)
    { res(iRes) = vec.elt(i);}
    return Rcpp::wrap(res);
  }
  /* support for wrap */
  template <typename Type, int SizeRows_, int SizeCols_, bool Orient_>
  SEXP wrap( STK::CArray<Type, SizeRows_, SizeCols_, Orient_> const& matrix)
  {
    enum
    {
      Rtype_ = STK::hidden::RcppTraits<Type>::Rtype_
    };

    Matrix<Rtype_> res(matrix.sizeRows(), matrix.sizeCols());
    for(int j=matrix.beginCols(), jRes=0; j< matrix.endCols(); j++, jRes++)
    {
      for(int i=matrix.beginRows(), iRes=0; i< matrix.endRows(); i++, iRes++)
      { res(iRes, jRes) = matrix.elt(i,j);}
    }
    return Rcpp::wrap(res);
  }
  /* support for wrap */
  template <typename Type, int SizeRows_, bool Orient_>
  SEXP wrap( STK::CArrayVector<Type, SizeRows_, Orient_ > const& vec)
  {
    enum
    {
      Rtype_ = STK::hidden::RcppTraits<Type>::Rtype_
    };

    Vector<Rtype_> res(vec.size());
    for(int i=vec.begin(), iRes=0; i< vec.end(); i++, iRes++)
    { res(iRes) = vec.elt(i);}
    return Rcpp::wrap(res);
  }
  template <typename Type, int SizeCols_, bool Orient_>
  SEXP wrap( STK::CArrayPoint<Type, SizeCols_, Orient_> const& vec)
  {
    enum
    {
      Rtype_ = STK::hidden::RcppTraits<Type>::Rtype_
    };
    Vector<Rtype_> res(vec.size());
    for(int i=vec.begin(), iRes=0; i< vec.end(); i++, iRes++)
    { res(iRes) = vec.elt(i);}
    return Rcpp::wrap(res);
  }
} // namespace Rcpp


#endif /* RCPP_WRAP_H */
