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
 * Project:  stkpp::StatDesc
 * Purpose:  Compute the confusion matrix
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Stat_ConfusionMatrix.h
 *  @brief This file contain the function confusionMatrix
 **/

#ifndef STK_STAT_CONFUSIONMATRIX_H
#define STK_STAT_CONFUSIONMATRIX_H

#include "STK_Stat_Factor.h"

namespace STK
{
namespace Stat
{

/** @ingroup StatDesc
 *  @brief Class computing the confusion matrix given two vectors
 *  of integer.
 **/
template < class TrueArray, class PredArray>
class ConfusionMatrix
{
  public:
    /** default constructor */
    ConfusionMatrix( TrueArray const& trueClass, PredArray const& predictClass)
                   : trueClass_(trueClass), predictClass_(predictClass)
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(TrueArray);
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(PredArray);
      if (trueClass_.range() != predictClass_.range())
      { STKRUNTIME_ERROR_NO_ARG(confusionMatrix(trueClass,predictClass),trueClass and predictClass are not of the same size);}
    }

    /** @return the confusion matrix obtained with the true classes and the
     *  predicted classes.
     **/
    CSquareXi operator()()
    {
      // compute nbLevels and create result
      int min = std::min(trueClass_.minElt(), predictClass_.minElt());
      int max = std::max(trueClass_.maxElt(), predictClass_.maxElt());
      CSquareXi res(Range(min, max, 0), 0);
      for (int i=trueClass_.begin(); i<trueClass_.end(); ++i)
      { res(trueClass_[i], predictClass_[i])++;}
      return res;
    }
  protected:
    TrueArray const& trueClass_;
    PredArray const& predictClass_;
};


/** @ingroup StatDesc
 *  @brief Utility function computing the confusion matrix given two vectors
 *  of integer.
 **/
template < class TrueArray, class PredArray>
CSquareXi confusionMatrix(TrueArray const& trueClass, PredArray const& predictClass)
{ return ConfusionMatrix<TrueArray, PredArray>(trueClass,predictClass)();}

}  // namespace Stat

}  // namespace STK

#endif /*STK_STAT_CONFUSIONMATRIX_H*/
