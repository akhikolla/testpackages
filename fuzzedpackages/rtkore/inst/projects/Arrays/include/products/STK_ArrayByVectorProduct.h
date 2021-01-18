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
 * Project:  stkpp::Arrays
 * created on: 30 déc. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ArrayByVectorProduct.h
 *  @brief In this file we implement the General Array by Vector product.
 **/


#ifndef STK_ARRAYBYVECTORPRODUCT_H
#define STK_ARRAYBYVECTORPRODUCT_H

namespace STK
{

namespace hidden
{

// point by array product
const int pointByArrayRowSize_ = 256;
const int pointByArrayColSize_ = 8;

typedef TRange<pointByArrayRowSize_> pointByArrayRowRange;
typedef TRange<pointByArrayColSize_> PointByArrayColRange;


/** @ingroup hidden
 *  this structure regroup all the methods using only pointers on the Type
 **/
template<typename Type>
struct MultImpl
{
  /** multiplication of two vectors */
  static Type vectorByVector(Type const* p_lhs, Type const* p_rhs)
  {
    Type sum = Type(0);
    for (int k=0; k< vectorSize; ++k) sum += p_lhs[k] * p_rhs[k];
    return(sum);
  }
  /** multiplication of two vectors */
  static Type PanelByVector(Type const* p_lhs, Type const* p_rhs)
  {
    Type sum = Type(0);
    for (int k=0; k< vectorSize; ++k) sum += p_lhs[k] * p_rhs[k];
    return(sum);
  }
};

/** @ingroup hidden
 *  Methods to use for C=AV with A a general matrix and V a vector.
 *  The structure bv contains only static methods and typedef and should normally
 *  not be used directly.
 **/
template<typename Lhs, typename Rhs, typename Result>
struct bv
{
  typedef typename Result::Type Type;
  /** Main method for Matrices by vector multiplication implementation
   *  res have been resized and initialized to zero outside
   *  this method.
   **/
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    for (int j= lhs.beginCols(); j<lhs.endCols(); ++j)
    {
      for (int i= lhs.beginRows(); i< lhs.endRows(); ++i)
      { res.elt(i) += lhs(i,j) * rhs[j];}
    }
    return;
  }
}; // struct bv

/** @ingroup hidden
 *  Methods to use for C=PB with P a point and B a matrix.
 *  The structure vb contains only static method and typedef and should normally
 *  not be used directly.
 **/
template<typename Lhs, typename Rhs, typename Result>
struct vb
{
  typedef typename Result::Type Type;
  /** Main method for point by Matrices multiplication implementation.
   *  res have been resized and initialized to zero outside
   *  this method.
   **/
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    for (int j= rhs.beginCols(); j<rhs.endCols(); ++j)
    {
      for (int i= rhs.beginRows(); i< rhs.endRows(); ++i)
      { res.elt(j) += lhs[i] * rhs(i,j);}
    }
  }
}; // struct pb

/** @ingroup hidden
 *  This structure regroup the products between a point and different kind of array
 **/
template<typename Lhs, typename Rhs, typename Result>
struct MultPointArray
{
  typedef typename Result::Type Type;
   /** Compute the product res = l*r with l a point (a row-vector) and r a matrix
   *   using the standard formula \f$ \mathrm{res}_j = \sum_{i=1}^n l_i r_{ij} \f$
   **/
  template<class SubLhs, class SubRhs>
  static void mult(SubLhs const& l, SubRhs const& r, Result& res)
  {
    // loop over the column of the right hand side
    for (int j=r.beginCols(); j< r.endCols(); ++j)
    {
      Type sum(0); // compute dot product
      for (int i=l.begin(); i<l.end(); ++i) { sum += l.elt(i) * r.elt(i,j); }
      res.elt(j) += sum;
    }
  }
  /** Compute the product res = l*r with l a point (a row-vector) and r a matrix
   *  using column decomposition.
   **/
  static void run( ExprBase<Lhs> const& l, ExprBase<Rhs> const& r, Result& res)
  {
    pointByArrayRowRange rowRange(r.beginRows(), pointByArrayRowSize_);
    for(;rowRange.end()<r.endRows(); rowRange += pointByArrayRowSize_)
    {
      PointByArrayColRange colRange(r.beginCols(), pointByArrayColSize_); // create fixed size col range
      for(; colRange.end() < r.endCols(); colRange += pointByArrayColSize_)
      {
        mult(l.sub(rowRange), r.sub(rowRange, colRange), res);
      }
      Range lastColRange(colRange.begin(), r.lastIdxCols(), 0);
      mult(l.sub(rowRange), r.sub(rowRange, lastColRange ), res);
    }
    // end
    Range lastRowRange(rowRange.begin(), r.lastIdxRows(), 0);
    PointByArrayColRange colRange(r.beginCols(), pointByArrayColSize_); // create fixed size col range
    for(; colRange.end() < r.endCols(); colRange += pointByArrayColSize_)
    { mult(l.sub(lastRowRange), r.sub(lastRowRange, colRange), res);}
    // last term
    Range lastColRange(colRange.begin(), r.lastIdxCols(), 0);
    mult(l.sub(lastRowRange), r.sub(lastRowRange, lastColRange ), res);
  }
};

} // namespace hidden



} // namespace STK

#endif /* STK_ARRAYBYVECTORPRODUCT_H */
