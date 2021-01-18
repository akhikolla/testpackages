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
 * Project: stkpp::STatistiK::StatDesc
 * Purpose: Perform the usual transformation on data set.
 * Author:  Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Stat_Transform.h
 *  @brief In this file we implement the main transformation on data set.
 **/

#ifndef STK_STAT_TRANSFORM_H
#define STK_STAT_TRANSFORM_H

#include <STKernel/include/STK_Misc.h>
#include "STK_Stat_Functors.h"

namespace STK
{
namespace Stat
{

/** @ingroup StatDesc
 *  Compute the mean by column of the variables in the container m and center it.
 *  @param m,mean data and means arrays
 **/
template < class Array, class RowVector>
void center( Array& m, RowVector& mean)
{
  typedef typename Array::Type Type;
  enum
  {
    sizeRows_ = Array::sizeRows_,
    sizeCols_ = Array::sizeCols_
  };
  m -= Const::Vector<Type, sizeRows_>(m.rows())
                * (mean = Stat::meanByCol(m));
}
template < class Array, class RowVector>
inline void centerByCol( Array& m, RowVector& mean)
{ center(m, mean);}

/** @ingroup StatDesc
 *  Compute the mean by row of the variables in the container V and center V.
 *  @param m,mean data and means arrays
 **/
template < class Array, class ColVector>
void centerByRow( Array& m, ColVector& mean)
{
  typedef typename Array::Type Type;
  enum
  {
    sizeRows_ = Array::sizeRows_,
    sizeCols_ = Array::sizeCols_
  };
  m -= (mean = Stat::meanByRow(m))
     * Const::Point<Type, sizeCols_>(m.cols());
}

/** @ingroup StatDesc
 *  Compute the weighted mean by column of the variables and center m.
 *  @param m,W,mean container of the data, weights and means
 **/
template < class Array, class Weights, class RowVector>
void center( Array& m, Weights const& W, RowVector& mean)
{
  typedef typename Array::Type Type;
  enum
  {
    sizeRows_ = Array::sizeRows_,
    sizeCols_ = Array::sizeCols_
  };
  m -= Const::Vector<Type, sizeRows_>(m.rows())
                 * (mean = Stat::meanByCol(m, W));
}
template < class Array, class Weights, class RowVector>
inline void centerByCol( Array& m, Weights const& W, RowVector& mean)
{ center(m, W, mean);}

/** @ingroup StatDesc
 *  Compute the weighted mean of the variables by rows and center m.
 *  @param m,W,mean container of the data, weights and means
 **/
template < class Array, class Weights, class ColVector>
void centerByRow( Array& m, Weights const& W, ColVector& mean)
{
  typedef typename Array::Type Type;
  enum
  {
    sizeRows_ = Array::sizeRows_,
    sizeCols_ = Array::sizeCols_
  };
  m -= (mean = Stat::meanByRow(m, W))
                 * Const::Point<Type, sizeCols_>(m.cols());
}

/** @ingroup StatDesc
 *  Compute the mean and the standard deviation by columns of the variable m
 *  and standardize it.
 *  @param m,mean,std arrays with the data, the means and the standard deviations
 *  @param unbiased @c false if the standard deviation as to be corrected,
 *  @c true otherwise
 **/
template < class Array,  class RowVector>
void standardize( Array& m, RowVector& mean, RowVector& std, bool unbiased = false)
{
  typedef typename Array::Type Type;
  enum
  {
    sizeRows_ = Array::sizeRows_,
    sizeCols_ = Array::sizeCols_
  };
  // center V
  centerByCol(m, mean);
  // compute variance with mean 0 as V is centered
  std = Stat::varianceWithFixedMeanByCol(m, Const::Vector<Type, sizeCols_>(m.cols())*Type(0), unbiased);
  for (int j= m.beginCols(); j< m.endCols(); j++)
  {
    Real dev= std[j];
    std[j] = (dev = Arithmetic<Real>::isFinite(dev) ? std::sqrt((double)dev) : 0.);
    if (dev) { m.col(j) /= dev;}
  }
}
template < class RowVector, class Array >
inline void standardizeByCol( Array& m, RowVector& mean, RowVector& std, bool unbiased = false)
{ standardize(m, mean, std, unbiased);}

/** @ingroup StatDesc
 *  Compute the mean and the standard deviation by rows of the variable m
 *  and standardize it.
 *  @param m,mean,std arrays with the data, the means and the standard deviations
 *  @param unbiased @c false if the standard deviation as to be corrected,
 *  @c true otherwise
 **/
template < class Array,  class ColVector>
void standardizeByRow( Array& m, ColVector& mean, ColVector& std, bool unbiased = false)
{
  typedef typename Array::Type Type;
  enum
  {
    sizeRows_ = Array::sizeRows_,
    sizeCols_ = Array::sizeCols_
  };
  // center V
  centerByRow(m, mean);
  // compute variance with mean 0 as V is centered
  std = Stat::varianceWithFixedMeanByRow(m, Const::Point<Type, sizeRows_>(m.rows())*Type(0), unbiased);
  for (int i= m.beginRows(); i< m.endRows(); i++)
  {
    Real dev = std[i];
    std[i] = (dev = Arithmetic<Real>::isFinite(dev) ? std::sqrt((double)dev) : 0.);
    if (dev) { m.row(i) /= dev;}
  }
}

/** @ingroup StatDesc
 *  Compute the weighted means and standard deviations by columns of the variable m
 *  and standardize it.
 *  @param m,W,mean,std containers with the Data, the weights, the means and standard deviations
 *  @param unbiased @c false if the standard deviation as to be corrected,
 *  @c true otherwise
 **/
template < class Array, class Weights, class RowVector>
void standardize( Array& m, Weights const& W, RowVector& mean, RowVector& std, bool unbiased = false)
{
  typedef typename Array::Type Type;
  enum
  {
    sizeRows_ = Array::sizeRows_,
    sizeCols_ = Array::sizeCols_
  };
  // center m
  centerByCol(m, W, mean);
  // compute variance with mean 0 as m is centered
  std = Stat::varianceWithFixedMeanByCol(m, W, Const::Vector<Type, sizeCols_>(m.cols())*Type(0), unbiased);
  // center
  for (int j= m.beginCols(); j< m.endCols(); j++)
  {
    // compute standard deviation
    Real dev = std[j];
    // take square root and save result
    std[j] = (dev = Arithmetic<Real>::isFinite(dev) ?  std::sqrt((double)dev) : 0.);
    // standardize data if necessary
    if (dev) { m.col(j) /= std[j];}
  }
}
template < class Array, class Weights, class RowVector>
inline void standardizeByCol( Array& m, Weights const& W, RowVector& mean, RowVector& std, bool unbiased = false)
{ standardize(m, W, mean, std, unbiased);}

/** @ingroup StatDesc
 *  Compute the weighted means and standard deviations by rows of the variable m
 *  and standardize it.
 *  @param m,W,mean,std containers with the Data, the weights, the means and standard deviations
 *  @param unbiased @c false if the standard deviation as to be corrected,
 *  @c true otherwise
 **/
template < class Array, class Weights, class ColVector>
void standardizeByRow( Array& m, Weights const& W, ColVector& mean, ColVector& std, bool unbiased = false)
{
  typedef typename Array::Type Type;
  enum
  {
    sizeRows_ = Array::sizeRows_,
    sizeCols_ = Array::sizeCols_
  };
  // center m
  centerByRow(m, W, mean);
  // compute variance with mean 0 as m is centered
  std = Stat::varianceWithFixedMeanByRow(m, W, Const::Point<Type, sizeRows_>(m.rows())*Type(0), unbiased);
  // center
  for (int i= m.beginRows(); i< m.endRows(); i++)
  {
    // compute standard deviation
    Real dev = std[i];
    // take square root and save result
    std[i] = (dev = Arithmetic<Real>::isFinite(dev) ?  std::sqrt((double)dev) : 0.);
    // standardize data if necessary
    if (dev) { m.row(i) /= std[i];}
  }
}

/** Add the means to the columns of the container m.
 *  @param m,mean the container with the data and the vector of the mean
 *  @param mean the Vector of the means
 **/
template < class RowVector, class Array >
void uncenter( Array&  m, RowVector const& mean)
{
  typedef typename Array::Type Type;
  enum
  {
    sizeRows_ = Array::sizeRows_,
    sizeCols_ = Array::sizeCols_
  };
  if (m.cols() != mean.range())
    STKRUNTIME_ERROR_NO_ARG(Stat::uncenter(m, mean),ranges are not the same);
  // apply uncenter
  m += Const::Vector<Type, sizeRows_>(m.rows())*mean;
}
template < class RowVector, class Array >
inline void uncenterByCol( Array&  m, RowVector const& mean)
{ uncenter(m, mean);}

/** Add the means to the rows of the container m.
 *  @param m,mean the container with the data and the vector of the mean
 *  @param mean the Vector of the means
 **/
template < class ColVector, class Array >
void uncenterByRow( Array&  m, ColVector const& mean)
{
  typedef typename Array::Type Type;
  enum
  {
    sizeRows_ = Array::sizeRows_,
    sizeCols_ = Array::sizeCols_
  };
  if (m.rows() != mean.range())
    STKRUNTIME_ERROR_NO_ARG(Stat::uncenterByRow(m, mean),ranges are not the same);
  // apply uncenter
  m += mean * Const::Point<Type, sizeCols_>(m.cols());
}
/** @ingroup StatDesc
 *  undo the standardization by columns of the standardized variable m.
 *  @param m,std the container with the data and the standard deviations
 **/
template < class Array, class RowVector>
void unstandardize( Array& m, RowVector const& std)
{
  typedef typename Array::Type Type;
  enum
  {
    sizeRows_ = Array::sizeRows_,
    sizeCols_ = Array::sizeCols_
  };
  if (m.cols() != std.range())
    STKRUNTIME_ERROR_NO_ARG(Error in Stat::unstandardize(m, std),ranges are not the sames);
  m = m.prod(Const::Vector<Type, sizeRows_>(m.rows())*std);
}
template < class Array, class RowVector>
inline void unstandardizeByCol( Array& m, RowVector const& std)
{ unstandardize(m, std);}

/** @ingroup StatDesc
 *  undo the standardization by rows of the standardized variable m.
 *  @param m,std the container with the data and the standard deviations
 **/
template < class Array, class ColVector>
void unstandardizeByRow( Array& m, ColVector const& std)
{
  typedef typename Array::Type Type;
  enum
  {
    sizeRows_ = Array::sizeRows_,
    sizeCols_ = Array::sizeCols_
  };
  if (m.rows() != std.range())
    STKRUNTIME_ERROR_NO_ARG(Error in Stat::unstandardizeByRow(m, std),ranges are not the sames);
  m = m.prod(std * Const::Point<Type, sizeCols_>(m.cols()));
}

/** @ingroup StatDesc
 *  undo the standardization by columns of the standardized variable m
 *  @param m,mean,std the containers with the Data, the means and the standard deviations
 **/
template < class Array, class RowVector>
void unstandardize( Array& m, RowVector const& mean, RowVector const& std)
{
  unstandardizeByCol(m, std); // unstandardize
  uncenterByCol(m, mean);     // uncenter
}
template < class Array, class RowVector>
inline void unstandardizeByCol( Array& m, RowVector const& mean, RowVector const& std)
{ unstandardize(m, mean, std);}

/** @ingroup StatDesc
 *  undo the standardization by rows of the standardized variable m
 *  @param m,mean,std the containers with the Data, the means and the standard deviations
 **/
template < class Array, class ColVector>
void unstandardizeByRow( Array& m, ColVector const& mean, ColVector const& std)
{
  unstandardizeByRow(m, std); // unstandardize
  uncenterByRow(m, mean);     // uncenter
}

}  // namespace Stat

}  // namespace STK

#endif /*STK_STAT_TRANSFORM_H*/
