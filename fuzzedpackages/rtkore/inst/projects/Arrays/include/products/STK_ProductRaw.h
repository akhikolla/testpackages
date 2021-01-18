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
 * created on: 25 déc. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ProductRaw.h
 *  @brief In this file we implement the raw static methods used by the products methods.
 **/

#ifndef STK_PRODUCTRAW_H
#define STK_PRODUCTRAW_H

namespace STK
{
/* size of the block and panels used in the product algorithm */
const int blockSize_ = 4;
const int panelSize_ = 64;
const int panelTotal= 256; // = 64 * 4

const int vectorSize = 256;

typedef TRange<panelSize_> panelRange;
typedef TRange<blockSize_> blockRange;

namespace hidden
{
/** @ingroup hidden
 *  This structure encapsulate the data allocated for a panel.
 **/
template<class Type>
struct Panel
{
  Type panel[blockSize_*panelSize_];
  inline Type const& operator[](int i) const { return panel[i];}
  inline Type& operator[](int i) { return panel[i];}
};

/** @ingroup hidden
 *  This structure encapsulate the data allocated for a block.
 **/
template<class Type>
struct Block
{
  Type block[blockSize_*blockSize_];
  inline Type const& operator[](int i) const { return block[i];}
  inline Type& operator[](int i) { return block[i];}
};

/** @ingroup hidden
 *  This structure encapsulate the data allocated for a block.
 **/
template<class Type>
struct RawVec
{
  Type vec[panelSize_];
  inline Type const& operator[](int i) const { return vec[i];}
  inline Type& operator[](int i) { return vec[i];}
};

/** @ingroup hidden
 *  This structure regroup the methods to used after block multiplication in
 *  order to perform the product of the remaining rows and columns.
 **/
template<typename Lhs, typename Rhs, typename Result>
struct MultCoefImpl
{
  typedef typename Result::Type Type;

  /** dot product. general by general*/
  static void dot( Lhs const& lhs, Rhs const& rhs, Result& res, int i, int j)
  {
    Range const dotRange = inf(lhs.rangeColsInRow(i), rhs.rangeRowsInCol(j));
    for (int k=dotRange.begin(); k< dotRange.end(); ++k)
      res.elt(i, j) += lhs.elt(i, k) * rhs.elt(k, j);
  }
  /** dot product. general by vector */
  static void dot( Lhs const& lhs, ITContainer<Rhs, Arrays::vector_> const& rhs, Result& res, int i)
  {
    Range const dotRange = inf(lhs.rangeColsInRow(i), rhs.range());
    for (int k=dotRange.begin(); k< dotRange.end(); ++k)
      res.elt(i) += lhs.elt(i, k) * rhs.elt(k);
  }
  /** dot product. point by general */
  static void dot( ITContainer<Lhs, Arrays::point_> const& lhs, Rhs const& rhs, Result& res, int j)
  {
    Range const dotRange = inf(rhs.rangeRowsInCol(j), lhs.range());
    for (int k=dotRange.begin(); k< dotRange.end(); ++k)
      res.elt(j) += lhs.elt(k) * rhs.elt(k, j);
  }

  static bool multDispatcher( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    // small sizes
    switch (lhs.sizeRows())
    {
      case 0: return true;
      case 1: mul1XX(lhs, rhs, res); return true;
      case 2: mul2XX(lhs, rhs, res); return true;
      case 3: mul3XX(lhs, rhs, res); return true;
      case 4: mul4XX(lhs, rhs, res); return true;
      case 5: mul5XX(lhs, rhs, res); return true;
      case 6: mul6XX(lhs, rhs, res); return true;
      case 7: mul7XX(lhs, rhs, res); return true;
      default: break;
    }
    switch (lhs.sizeCols())
    {
      case 0: return true; break;
      case 1: mulX1X(lhs, rhs, res); return true;
      case 2: mulX2X(lhs, rhs, res); return true;
      case 3: mulX3X(lhs, rhs, res); return true;
      case 4: mulX4X(lhs, rhs, res); return true;
      case 5: mulX5X(lhs, rhs, res); return true;
      case 6: mulX6X(lhs, rhs, res); return true;
      case 7: mulX7X(lhs, rhs, res); return true;
      default: break;
    }
    switch (rhs.sizeCols())
    {
      case 0: return true;
      case 1: mulXX1(lhs, rhs, res); return true;
      case 2: mulXX2(lhs, rhs, res); return true;
      case 3: mulXX3(lhs, rhs, res); return true;
      case 4: mulXX4(lhs, rhs, res); return true;
      case 5: mulXX5(lhs, rhs, res); return true;
      case 6: mulXX6(lhs, rhs, res); return true;
      case 7: mulXX7(lhs, rhs, res); return true;
      default: break;
    }
    return false;
  }
  /** multiplication with Lhs having 1 row */
  static void mul1XX( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const lhsRow = lhs.beginRows();
    for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      for (int k=rhs.beginRows(); k< rhs.endRows(); ++k)
        res.elt(lhsRow, j) += lhs.elt(lhsRow, k) * rhs.elt(k, j);
  }
  /** multiplication with Lhs having 2 rows */
  static void mul2XX( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const i = lhs.beginRows();
    for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      for (int k=rhs.beginRows(); k< rhs.endRows(); ++k)
       {
         res.elt(i  , j) += lhs.elt(i  , k) * rhs.elt(k, j);
         res.elt(i+1, j) += lhs.elt(i+1, k) * rhs.elt(k, j);
       }
  }
  /** multiplication with Lhs having 3 rows */
  static void mul3XX( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const i = lhs.beginRows();
    for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      for (int k=rhs.beginRows(); k< rhs.endRows(); ++k)
      {
        res.elt(i  , j) += lhs.elt(i  , k) * rhs.elt(k, j);
        res.elt(i+1, j) += lhs.elt(i+1, k) * rhs.elt(k, j);
        res.elt(i+2, j) += lhs.elt(i+2, k) * rhs.elt(k, j);
      }
  }
  /** multiplication with Lhs having 4 rows */
  static void mul4XX( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const i = lhs.beginRows();
    for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      for (int k=rhs.beginRows(); k< rhs.endRows(); ++k)
      {
        res.elt(i  , j) += lhs.elt(i  , k) * rhs.elt(k, j);
        res.elt(i+1, j) += lhs.elt(i+1, k) * rhs.elt(k, j);
        res.elt(i+2, j) += lhs.elt(i+2, k) * rhs.elt(k, j);
        res.elt(i+3, j) += lhs.elt(i+3, k) * rhs.elt(k, j);
      }
  }
  /** multiplication with Lhs having 5 rows */
  static void mul5XX( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const i = lhs.beginRows();
    for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      for (int k=rhs.beginRows(); k< rhs.endRows(); ++k)
      {
        res.elt(i  , j) += lhs.elt(i  , k) * rhs.elt(k, j);
        res.elt(i+1, j) += lhs.elt(i+1, k) * rhs.elt(k, j);
        res.elt(i+2, j) += lhs.elt(i+2, k) * rhs.elt(k, j);
        res.elt(i+3, j) += lhs.elt(i+3, k) * rhs.elt(k, j);
        res.elt(i+4, j) += lhs.elt(i+4, k) * rhs.elt(k, j);
      }
  }
  /** multiplication with Lhs having 6 rows */
  static void mul6XX( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const i = lhs.beginRows();
    for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      for (int k=rhs.beginRows(); k< rhs.endRows(); ++k)
      {
        res.elt(i  , j) += lhs.elt(i  , k) * rhs.elt(k, j);
        res.elt(i+1, j) += lhs.elt(i+1, k) * rhs.elt(k, j);
        res.elt(i+2, j) += lhs.elt(i+2, k) * rhs.elt(k, j);
        res.elt(i+3, j) += lhs.elt(i+3, k) * rhs.elt(k, j);
        res.elt(i+4, j) += lhs.elt(i+4, k) * rhs.elt(k, j);
        res.elt(i+5, j) += lhs.elt(i+5, k) * rhs.elt(k, j);
      }
  }
  /** multiplication with Lhs having 7 rows */
  static void mul7XX( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const i = lhs.beginRows();
    for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      for (int k=rhs.beginRows(); k< rhs.endRows(); ++k)
      {
        res.elt(i  , j) += lhs.elt(i  , k) * rhs.elt(k, j);
        res.elt(i+1, j) += lhs.elt(i+1, k) * rhs.elt(k, j);
        res.elt(i+2, j) += lhs.elt(i+2, k) * rhs.elt(k, j);
        res.elt(i+3, j) += lhs.elt(i+3, k) * rhs.elt(k, j);
        res.elt(i+4, j) += lhs.elt(i+4, k) * rhs.elt(k, j);
        res.elt(i+5, j) += lhs.elt(i+5, k) * rhs.elt(k, j);
        res.elt(i+6, j) += lhs.elt(i+6, k) * rhs.elt(k, j);
      }
  }
  /** multiplication with Lhs having 1 column and Rhs having 1 row */
  static void mulX1X( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const k = lhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
        res.elt(i, j) += lhs.elt(i, k) * rhs.elt(k, j);
  }
  /** multiplication with Lhs having 2 columns and Rhs having 2 rows */
  static void mulX2X( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const k = lhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      {
        res.elt(i, j) += lhs.elt(i, k)   * rhs.elt(k, j);
        res.elt(i, j) += lhs.elt(i, k+1) * rhs.elt(k+1, j);
      }
  }
  /** multiplication with Lhs having 3 columns and Rhs having 3 rows */
  static void mulX3X( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const k = lhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      {
        res.elt(i, j) += lhs.elt(i, k)   * rhs.elt(k, j);
        res.elt(i, j) += lhs.elt(i, k+1) * rhs.elt(k+1, j);
        res.elt(i, j) += lhs.elt(i, k+2) * rhs.elt(k+2, j);
      }
  }
  /** multiplication with Lhs having 4 columns and Rhs having 4 rows */
  static void mulX4X( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const k = lhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      {
        res.elt(i, j) += lhs.elt(i, k)   * rhs.elt(k, j);
        res.elt(i, j) += lhs.elt(i, k+1) * rhs.elt(k+1, j);
        res.elt(i, j) += lhs.elt(i, k+2) * rhs.elt(k+2, j);
        res.elt(i, j) += lhs.elt(i, k+3) * rhs.elt(k+3, j);
      }
  }
  /** multiplication with Lhs having 5 columns and Rhs having 5 rows */
  static void mulX5X( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const k = lhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      {
        res.elt(i, j) += lhs.elt(i, k)   * rhs.elt(k, j);
        res.elt(i, j) += lhs.elt(i, k+1) * rhs.elt(k+1, j);
        res.elt(i, j) += lhs.elt(i, k+2) * rhs.elt(k+2, j);
        res.elt(i, j) += lhs.elt(i, k+3) * rhs.elt(k+3, j);
        res.elt(i, j) += lhs.elt(i, k+4) * rhs.elt(k+4, j);
      }
  }
  /** multiplication with Lhs having 6 columns and Rhs having 6 rows */
  static void mulX6X( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const k = lhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      {
        res.elt(i, j) += lhs.elt(i, k)   * rhs.elt(k, j);
        res.elt(i, j) += lhs.elt(i, k+1) * rhs.elt(k+1, j);
        res.elt(i, j) += lhs.elt(i, k+2) * rhs.elt(k+2, j);
        res.elt(i, j) += lhs.elt(i, k+3) * rhs.elt(k+3, j);
        res.elt(i, j) += lhs.elt(i, k+4) * rhs.elt(k+4, j);
        res.elt(i, j) += lhs.elt(i, k+5) * rhs.elt(k+5, j);
      }
  }
  /** multiplication with Lhs having 7 columns and Rhs having 7 rows */
  static void mulX7X( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const k = lhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      {
        res.elt(i, j) += lhs.elt(i, k)   * rhs.elt(k, j);
        res.elt(i, j) += lhs.elt(i, k+1) * rhs.elt(k+1, j);
        res.elt(i, j) += lhs.elt(i, k+2) * rhs.elt(k+2, j);
        res.elt(i, j) += lhs.elt(i, k+3) * rhs.elt(k+3, j);
        res.elt(i, j) += lhs.elt(i, k+4) * rhs.elt(k+4, j);
        res.elt(i, j) += lhs.elt(i, k+5) * rhs.elt(k+5, j);
        res.elt(i, j) += lhs.elt(i, k+6) * rhs.elt(k+6, j);
       }
  }
  /** multiplication with Rhs having 1 column */
  static void mulXX1( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const j = rhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int k=lhs.beginCols(); k< lhs.endCols(); ++k)
        res.elt(i, j) += lhs.elt(i, k) * rhs.elt(k, j);
  }
  /** multiplication with Rhs having 2 columns */
  static void mulXX2( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const j = rhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int k=lhs.beginCols(); k< lhs.endCols(); ++k)
      {
        res.elt(i, j  ) += lhs.elt(i, k) * rhs.elt(k, j  );
        res.elt(i, j+1) += lhs.elt(i, k) * rhs.elt(k, j+1);
      }
  }
  /** multiplication with Rhs having 3 columns */
  static void mulXX3( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const j = rhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int k=lhs.beginCols(); k< lhs.endCols(); ++k)
      {
        res.elt(i, j  ) += lhs.elt(i, k) * rhs.elt(k, j  );
        res.elt(i, j+1) += lhs.elt(i, k) * rhs.elt(k, j+1);
        res.elt(i, j+2) += lhs.elt(i, k) * rhs.elt(k, j+2);
      }
  }
  /** multiplication with Rhs having 4 columns */
  static void mulXX4( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const j = rhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int k=lhs.beginCols(); k< lhs.endCols(); ++k)
      {
        res.elt(i, j  ) += lhs.elt(i, k) * rhs.elt(k, j  );
        res.elt(i, j+1) += lhs.elt(i, k) * rhs.elt(k, j+1);
        res.elt(i, j+2) += lhs.elt(i, k) * rhs.elt(k, j+2);
        res.elt(i, j+3) += lhs.elt(i, k) * rhs.elt(k, j+3);
      }
  }
  /** multiplication with Rhs having 5 columns */
  static void mulXX5( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const j = rhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int k=lhs.beginCols(); k< lhs.endCols(); ++k)
      {
        res.elt(i, j  ) += lhs.elt(i, k) * rhs.elt(k, j  );
        res.elt(i, j+1) += lhs.elt(i, k) * rhs.elt(k, j+1);
        res.elt(i, j+2) += lhs.elt(i, k) * rhs.elt(k, j+2);
        res.elt(i, j+3) += lhs.elt(i, k) * rhs.elt(k, j+3);
        res.elt(i, j+4) += lhs.elt(i, k) * rhs.elt(k, j+4);
      }
  }
  /** multiplication with Rhs having 6 columns */
  static void mulXX6( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const j = rhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int k=lhs.beginCols(); k< lhs.endCols(); ++k)
      {
        res.elt(i, j  ) += lhs.elt(i, k) * rhs.elt(k, j  );
        res.elt(i, j+1) += lhs.elt(i, k) * rhs.elt(k, j+1);
        res.elt(i, j+2) += lhs.elt(i, k) * rhs.elt(k, j+2);
        res.elt(i, j+3) += lhs.elt(i, k) * rhs.elt(k, j+3);
        res.elt(i, j+4) += lhs.elt(i, k) * rhs.elt(k, j+4);
        res.elt(i, j+5) += lhs.elt(i, k) * rhs.elt(k, j+5);
      }
  }
  /** multiplication with Rhs having 7 columns */
  static void mulXX7( Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int const j = rhs.beginCols();
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int k=lhs.beginCols(); k< lhs.endCols(); ++k)
      {
        res.elt(i, j  ) += lhs.elt(i, k) * rhs.elt(k, j  );
        res.elt(i, j+1) += lhs.elt(i, k) * rhs.elt(k, j+1);
        res.elt(i, j+2) += lhs.elt(i, k) * rhs.elt(k, j+2);
        res.elt(i, j+3) += lhs.elt(i, k) * rhs.elt(k, j+3);
        res.elt(i, j+4) += lhs.elt(i, k) * rhs.elt(k, j+4);
        res.elt(i, j+5) += lhs.elt(i, k) * rhs.elt(k, j+5);
        res.elt(i, j+6) += lhs.elt(i, k) * rhs.elt(k, j+6);
      }
  }
  /** multiplication with one outer rows/columns */
  static void mult1Outer( Lhs const& lhs, Rhs const& rhs, Result& res, int k)
  {
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
        res.elt(i,j) += lhs.elt(i, k) * rhs.elt(k, j);
  }
  /** multiplication with two outer rows/columns */
  static void mult2Outer( Lhs const& lhs, Rhs const& rhs, Result& res, int k)
  {
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
        res.elt(i,j) += lhs.elt(i, k  ) * rhs.elt(k,   j)
                      + lhs.elt(i, k+1) * rhs.elt(k+1, j);
  }
  /** multiplication with three outer rows/columns */
  static void mult3Outer( Lhs const& lhs, Rhs const& rhs, Result& res, int k)
  {
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
        res.elt(i,j) += lhs.elt(i, k  ) * rhs.elt(k, j)
                      + lhs.elt(i, k+1) * rhs.elt(k+1, j)
                      + lhs.elt(i, k+2) * rhs.elt(k+2, j);
  }
};

/** @ingroup hidden
 *  This structure regroup the methods to used for copying part of an array in a Block or a Panel
 **/
template<typename Array, typename Type = typename Array::Type>
struct CopySubArrayImpl
{
  /**  copy a part of an Array in a Block */
  static void arrayToBlock( Array const& m, Block<Type>& block, int iRow, int jCol)
  {
     block[0]  = m.elt(iRow  , jCol);
     block[1]  = m.elt(iRow  , jCol+1);
     block[2]  = m.elt(iRow  , jCol+2);
     block[3]  = m.elt(iRow  , jCol+3);
     block[4]  = m.elt(iRow+1, jCol);
     block[5]  = m.elt(iRow+1, jCol+1);
     block[6]  = m.elt(iRow+1, jCol+2);
     block[7]  = m.elt(iRow+1, jCol+3);
     block[8]  = m.elt(iRow+2, jCol);
     block[9]  = m.elt(iRow+2, jCol+1);
     block[10] = m.elt(iRow+2, jCol+2);
     block[11] = m.elt(iRow+2, jCol+3);
     block[12] = m.elt(iRow+3, jCol);
     block[13] = m.elt(iRow+3, jCol+1);
     block[14] = m.elt(iRow+3, jCol+2);
     block[15] = m.elt(iRow+3, jCol+3);
  }
  /** copy a part of an Array in a Block with block columns size given */
  static void arrayToBlock( Array const& m, Block<Type>& block, int iRow, int jCol, int bSize)
  {
    for (int i=0; i<bSize; ++i)
    {
      block[i*blockSize_  ] = m.elt(iRow+i, jCol);
      block[i*blockSize_+1] = m.elt(iRow+i, jCol+1);
      block[i*blockSize_+2] = m.elt(iRow+i, jCol+2);
      block[i*blockSize_+3] = m.elt(iRow+i, jCol+3);
    }
  }
  /** copy a part of an array in a Panel */
  static void arrayToPanel( Array const& m, Panel<Type>& panel, int iRow, int jCol)
  {
    for (int j=0; j<panelSize_; ++j)
    {
      panel[j*blockSize_  ] = m.elt(iRow  , jCol+j);
      panel[j*blockSize_+1] = m.elt(iRow+1, jCol+j);
      panel[j*blockSize_+2] = m.elt(iRow+2, jCol+j);
      panel[j*blockSize_+3] = m.elt(iRow+3, jCol+j);
    }
  }
  /** copy a part of an array in a panel with Panel columns size given */
  static void arrayToPanel( Array const& rhs, Panel<Type>& panel, int iRow, int jCol, int pSize)
  {
    for (int j=0; j<pSize; ++j)
    {
      panel[j*blockSize_  ] = rhs.elt(iRow  , jCol+j);
      panel[j*blockSize_+1] = rhs.elt(iRow+1, jCol+j);
      panel[j*blockSize_+2] = rhs.elt(iRow+2, jCol+j);
      panel[j*blockSize_+3] = rhs.elt(iRow+3, jCol+j);
    }
  }
  /** default dimensions */
  static void arrayToBlockByCol( Array const& m, Block<Type>& block, int iRow, int jCol)
  {
    block[0]  = m.elt(iRow  , jCol);
    block[1]  = m.elt(iRow+1, jCol);
    block[2]  = m.elt(iRow+2, jCol);
    block[3]  = m.elt(iRow+3, jCol);
    block[4]  = m.elt(iRow  , jCol+1);
    block[5]  = m.elt(iRow+1, jCol+1);
    block[6]  = m.elt(iRow+2, jCol+1);
    block[7]  = m.elt(iRow+3, jCol+1);
    block[8]  = m.elt(iRow  , jCol+2);
    block[9]  = m.elt(iRow+1, jCol+2);
    block[10] = m.elt(iRow+2, jCol+2);
    block[11] = m.elt(iRow+3, jCol+2);
    block[12] = m.elt(iRow  , jCol+3);
    block[13] = m.elt(iRow+1, jCol+3);
    block[14] = m.elt(iRow+2, jCol+3);
    block[15] = m.elt(iRow+3, jCol+3);
  }
  /** with block size given */
  static void arrayToBlockByCol( Array const& m, Block<Type>& block, int iRow, int jCol, int bSize)
  {
    for (int j=0; j<bSize; ++j)
    {
      block[j*blockSize_  ] = m.elt(iRow,   jCol+j);
      block[j*blockSize_+1] = m.elt(iRow+1, jCol+j);
      block[j*blockSize_+2] = m.elt(iRow+2, jCol+j);
      block[j*blockSize_+3] = m.elt(iRow+3, jCol+j);
    }
  }
  /** default dimensions */
  static void arrayToPanelByCol( Array const& m, Panel<Type>& panel, int iRow, int kPos)
  {
    for (int i=0; i<panelSize_; ++i)
    {
      panel[i*blockSize_  ] = m.elt(iRow+i,kPos);
      panel[i*blockSize_+1] = m.elt(iRow+i,kPos+1);
      panel[i*blockSize_+2] = m.elt(iRow+i,kPos+2);
      panel[i*blockSize_+3] = m.elt(iRow+i,kPos+3);
    }
  }
  /** with panel size dimension given */
  static void arrayToPanelByCol( Array const& m, Panel<Type>& panel, int iRow, int kPos, int pSize)
  {
    for (int i=0; i<pSize; ++i)
    {
      panel[i*blockSize_  ] = m.elt(iRow+i, kPos);
      panel[i*blockSize_+1] = m.elt(iRow+i, kPos+1);
      panel[i*blockSize_+2] = m.elt(iRow+i, kPos+2);
      panel[i*blockSize_+3] = m.elt(iRow+i, kPos+3);
    }
  }

};

} // namespace hidden

} // namespace STK

#endif /* STK_PRODUCTRAW_H */
