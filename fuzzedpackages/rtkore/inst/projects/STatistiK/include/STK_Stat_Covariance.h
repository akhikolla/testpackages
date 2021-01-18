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
 * Project:  stkpp::StatistiK::StatDesc
 * Purpose:  Compute elementary statistics for two variables.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Stat_Covariance.h
 *  @brief This file contains the methods computing the covariance of an array
 **/

#ifndef STK_STAT_COVARIANCE_H
#define STK_STAT_COVARIANCE_H

#include "STK_Stat_Functors.h"

namespace STK
{
namespace Stat
{

/** @ingroup StatDesc
 *  Compute the covariance of the variables X and Y.
 *  \f[ \hat{cov}(X,Y) = \frac{1}{n}
 *                  \sum_{i=1}^n (X_i - \bar{X})(Y_i - \bar{Y}).
 *  \f]
 *  @param X,Y variables
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template<class XArray, class YArray>
Real covariance( ExprBase<XArray> const& X
               , ExprBase<YArray> const& Y
               , bool unbiased = false
               )
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(XArray);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(YArray);
  typedef typename hidden::Traits<XArray>::Type XType;
  typedef typename hidden::Traits<YArray>::Type YType;
  typedef typename hidden::Promote<XType, YType>::result_type Type;

#ifdef STK_BOUNDS_CHECK
  if (X.range() != Y.range())
      STKRUNTIME_ERROR_NO_ARG(Error in Stat::covariance(X,Y,unbiased),ranges are not the sames);
#endif
  Type xMean = mean(X.asDerived()), yMean = mean(Y.asDerived());
  int nobs = X.size();
  Type xsum  = 0.0, ysum = 0.0, cov  = 0.0, xdev, ydev;
  for (int i=X.begin(); i<X.end(); i++)
  {
    xsum += (xdev = X[i] - xMean);
    ysum += (ydev = Y[i] - yMean);
    cov += (xdev*ydev);
  }
  // compute the covariance with corrected bias
  return (unbiased) ?
    (nobs > 1) ? ((cov - (xsum*ysum)/Type(nobs))/Type(nobs -1)) : 0.

    :
    (nobs > 0) ? (cov - (xsum*ysum)/Type(nobs))/Type(nobs) : 0.;
}


/** @ingroup StatDesc
 *  Compute the covariance of the variables X and Y
 *  \f[ \hat{cov}(X,Y) = \frac{1}{n}
 *                  \sum_{i=1}^n (X_i - \hat{X})(Y_i - \hat{Y}).
 *  \f]
 *  @param X,Y variables
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template<class XArray, class YArray>
Real covarianceSafe( ExprBase<XArray> const& X
                   , ExprBase<YArray> const& Y
                   , bool unbiased = false
                   )
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(XArray);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(YArray);
  typedef typename hidden::Traits<XArray>::Type XType;
  typedef typename hidden::Traits<YArray>::Type YType;
  typedef typename hidden::Promote<XType, YType>::result_type Type;

  #ifdef STK_BOUNDS_CHECK
    if (X.range() != Y.range())
        STKRUNTIME_ERROR_NO_ARG(Error in Stat::covarianceSafe(X,Y,unbiased),ranges are not the sames);
  #endif

  int nobs = X.size();
  Type xMean = mean(X.asDerived()), yMean = mean(Y.asDerived());
  Type xsum  = 0.0, ysum = 0.0, cov  = 0.0, xdev, ydev;
  for (int i=X.begin(); i<X.end(); i++)
  {
    if (Arithmetic<Type>::isFinite(X[i]) && Arithmetic<Type>::isFinite(Y[i]))
    {
      xsum += (xdev = X[i] - xMean);
      ysum += (ydev = Y[i] - yMean);
      cov += (xdev*ydev);
    }
    else nobs--;
  }
  // compute the covariance with corrected bias
  return (unbiased) ?
    (nobs > 1) ? ((cov - (xsum*ysum)/Type(nobs))/Type(nobs -1)) : 0.

    :
    (nobs > 0) ? (cov - (xsum*ysum)/Type(nobs))/Type(nobs) : 0.;
}

/** @ingroup StatDesc
 *
 *  The covariance between X and Y will be computed using the formula
 *  @see Mark Galassi, Jim Davies, James Theiler,
 *  Brian Gough, Gerard Jungman, Michael Booth, and Fabrice Rossi.
 *  http://www.gnu.org/software/gsl/manual GNU Scientific Library -
 *   Reference manual, Version 1.15, 2011.
 *  http://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.html
 *  Sec. 21.7 Weighted Samples
 *
 *  \f$ cov(X,Y)=\frac{\sum_{i=1}^{N}w_i}{\left(\sum_{i=1}^{N}w_i\right)^2-\sum_{i=1}^{N}w_i^2}
 *  \sum_{i=1}^N w_i \left(  x_{i}-\bar{x} \right)  \left( y_{i}-\bar{y} \right). \f$
 *
 *  @param X,Y,W variables and weights
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template<class XArray, class YArray, class Weights>
Real covariance( ExprBase<XArray> const& X
               , ExprBase<YArray> const& Y
               , ExprBase<Weights> const& W
               , bool unbiased = false
               )
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(XArray);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(YArray);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
  typedef typename hidden::Traits<XArray>::Type XType;
  typedef typename hidden::Traits<YArray>::Type YType;
  typedef typename hidden::Promote<XType, YType>::result_type Type;

#ifdef STK_BOUNDS_CHECK
  if (X.range() != Y.range())
      STKRUNTIME_ERROR_NO_ARG(Error in Stat::covarianceWithFixedMeanSafe(X,Y,W,xMu,yMu,unbiased),ranges are not the sames);
#endif

  Type xMean = mean(X.asDerived(),W.asDerived()), yMean = mean(Y.asDerived(),W.asDerived());
  Real xsum  = 0.0, ysum = 0.0, xdev, ydev, sumWeights = 0.0, sum2Weights = 0.0, cov = 0.0;
  for (int i=X.begin(); i<X.end(); i++)
  {
    xsum += (xdev = X[i] - xMean);
    ysum += (ydev = Y[i] - yMean);
    Real Wi = std::abs((Real)W[i]);
    cov         += Wi * (xdev*ydev);
    sumWeights  += Wi;
    sum2Weights += Wi * Wi;
  }
  // compute the variance
  return (unbiased) ?
       (sumWeights*sumWeights > sum2Weights) ?
         (cov - xsum*ysum)/(sumWeights - sum2Weights/sumWeights) : 0.
//         (cov - xsum*ysum/sumWeights)/(sumWeights - sum2Weights/sumWeights) : 0.
       : (sumWeights) ? (cov - xsum*ysum)/(sumWeights) : 0.;
}

/** @ingroup StatDesc
 *
 *  The covariance between X and Y will be computed using the formula
 *  @see Mark Galassi, Jim Davies, James Theiler,
 *  Brian Gough, Gerard Jungman, Michael Booth, and Fabrice Rossi.
 *  http://www.gnu.org/software/gsl/manual GNU Scientific Library -
 *   Reference manual, Version 1.15, 2011.
 *  http://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.html
 *  Sec. 21.7 Weighted Samples
 *
 *  \f$ cov(X,Y)=\frac{\sum_{i=1}^{N}w_i}{\left(\sum_{i=1}^{N}w_i\right)^2-\sum_{i=1}^{N}w_i^2}
 *  \sum_{i=1}^N w_i \left(  x_{i}-\bar{x} \right)  \left( y_{i}-\bar{y} \right). \f$
 *
 *  @param X,Y,W variables and weights
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template<class XArray, class YArray, class Weights>
Real covarianceSafe( ExprBase<XArray> const& X
                   , ExprBase<YArray> const& Y
                   , ExprBase<Weights> const& W
                   , bool unbiased = false
                   )
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(XArray);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(YArray);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
  typedef typename hidden::Traits<XArray>::Type XType;
  typedef typename hidden::Traits<YArray>::Type YType;
  typedef typename hidden::Promote<XType, YType>::result_type Type;

#ifdef STK_BOUNDS_CHECK
  if (X.range() != Y.range())
      STKRUNTIME_ERROR_NO_ARG(Error in Stat::covarianceWithFixedMeanSafe(X,Y,W,xMu,yMu,unbiased),ranges are not the sames);
#endif

  Type xMean = mean(X.asDerived(),W.asDerived()), yMean = mean(Y.asDerived(),W.asDerived());
  // compute covariance
  Type xsum  = 0.0, ysum = 0.0, xdev, ydev, sumWeights = 0.0, sum2Weights = 0.0, cov = 0.0;
  for (int i=X.begin(); i<X.end(); i++)
  {
    if  (Arithmetic<Type>::isFinite(X[i])
      && Arithmetic<Type>::isFinite(Y[i])
      && Arithmetic<Type>::isFinite(W[i])
      )
    {
      xsum += (xdev = X[i] - xMean); // deviation from the mean
      ysum += (ydev = Y[i] - yMean); // deviation from the mean
      Type Wi = std::abs((Real)W[i]);
      cov += Wi * (xdev*ydev);         // cross product
      sumWeights  += Wi;           // sum absolute weights
      sum2Weights += Wi * Wi;     // sum squared weights
    }
  }
  // compute the variance
  return (unbiased) ?
       (sumWeights*sumWeights > sum2Weights) ?
         (cov - xsum*ysum)/(sumWeights - sum2Weights/sumWeights) : 0.
//         (cov - xsum*ysum/sumWeights)/(sumWeights - sum2Weights/sumWeights) : 0.
       : (sumWeights) ? (cov - xsum*ysum)/(sumWeights) : 0.;
}

/** @ingroup StatDesc
 *  Compute the covariance of the variables X and Y with fixed
 *  means.
 *  \f[ \hat{cov}(X,Y) = \frac{1}{n}
 *                  \sum_{i=1}^n (X_i - \mu_X)(Y_i - \mu_Y).
 *  \f]
 *  @param X,Y variables
 *  @param xMu,yMu means of the variables
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template<class XArray, class YArray>
Real covarianceWithFixedMean( ExprBase<XArray> const& X
                            , ExprBase<YArray> const& Y
                            , typename hidden::Traits<XArray>::Type const& xMu
                            , typename hidden::Traits<YArray>::Type const& yMu
                            , bool unbiased = false
                            )
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(XArray);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(YArray);
  typedef typename hidden::Traits<XArray>::Type XType;
  typedef typename hidden::Traits<YArray>::Type YType;
  typedef typename hidden::Promote<XType, YType>::result_type Type;

#ifdef STK_BOUNDS_CHECK
  if (X.range() != Y.range())
      STKRUNTIME_ERROR_NO_ARG(Error in Stat::covarianceWithFixedMean(X,Y,xMu,yMu,unbiased),ranges are not the sames);
#endif

  int nobs = X.size();
  Type xsum  = 0.0, ysum = 0.0, cov  = 0.0, xdev, ydev;
  for (int i=X.begin(); i<X.end(); i++)
  {
    xsum += (xdev = X[i] - xMu);
    ysum += (ydev = Y[i] - yMu);
    cov += (xdev*ydev);
  }
  // compute the covariance with corrected bias
  return (unbiased) ?
    (nobs > 1) ? ((cov - (xsum*ysum)/Type(nobs))/Type(nobs -1))
               : Arithmetic<Real>::infinity()

    :
    (nobs > 0) ? (cov - (xsum*ysum)/Type(nobs))/Type(nobs)
               : Arithmetic<Real>::infinity();
}

/** @ingroup StatDesc
 *  Compute the covariance of the variables X and Y with fixed
 *  means.
 *  \f[ \hat{cov}(X,Y) = \frac{1}{n}
 *                  \sum_{i=1}^n (X_i - \mu_X)(Y_i - \mu_Y).
 *  \f]
 *  @param X,Y variables
 *  @param xMu,yMu means of the variables
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template<class XArray, class YArray>
Real covarianceWithFixedMeanSafe( ExprBase<XArray> const& X
                                , ExprBase<YArray> const& Y
                                , typename hidden::Traits<XArray>::Type const& xMu
                                , typename hidden::Traits<YArray>::Type const& yMu
                                , bool unbiased = false
                                )
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(XArray);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(YArray);
  typedef typename hidden::Traits<XArray>::Type XType;
  typedef typename hidden::Traits<YArray>::Type YType;
  typedef typename hidden::Promote<XType, YType>::result_type Type;

  #ifdef STK_BOUNDS_CHECK
    if (X.range() != Y.range())
        STKRUNTIME_ERROR_NO_ARG(Error in Stat::covarianceWithFixedMeanSafe(X,Y,xMu,yMu,unbiased),ranges are not the sames);
  #endif

  int nobs = X.size();
  Type xsum  = 0.0, ysum = 0.0, cov  = 0.0, xdev, ydev;
  for (int i=X.begin(); i<X.end(); i++)
  {
    if (Arithmetic<Type>::isFinite(X[i]) && Arithmetic<Type>::isFinite(Y[i]))
    {
      xsum += (xdev = X[i] - xMu); // deviation from the mean
      ysum += (ydev = Y[i] - yMu); // deviation from the mean
      cov += (xdev*ydev);         // squared value
    }
    else nobs--;
  }
  // compute the covariance with corrected bias
  return (unbiased) ?
    (nobs > 1) ? ((cov - (xsum*ysum)/Type(nobs))/Type(nobs -1))
               : Arithmetic<Real>::infinity()

    :
    (nobs > 0) ? (cov - (xsum*ysum)/Type(nobs))/Type(nobs)
               : Arithmetic<Real>::infinity();
}

/** @ingroup StatDesc
 *
 *  The covariance between X and Y will be computed using the formula
 *  @see Mark Galassi, Jim Davies, James Theiler,
 *  Brian Gough, Gerard Jungman, Michael Booth, and Fabrice Rossi.
 *  http://www.gnu.org/software/gsl/manual GNU Scientific Library -
 *   Reference manual, Version 1.15, 2011.
 *  http://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.html
 *  Sec. 21.7 Weighted Samples
 *
 *  \f$ cov(X,Y)=\frac{\sum_{i=1}^{N}w_i}{\left(\sum_{i=1}^{N}w_i\right)^2-\sum_{i=1}^{N}w_i^2}
 *  \sum_{i=1}^N w_i \left(  x_{i}-\bar{x} \right)  \left( y_{i}-\bar{y} \right). \f$
 *
 *  @param X,Y,W variables and weights
 *  @param xMean,yMean means of the varaibles
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template<class XArray, class YArray, class Weights>
Real covarianceWithFixedMean( ExprBase<XArray> const& X
                            , ExprBase<YArray> const& Y
                            , ExprBase<Weights> const& W
                            , typename hidden::Traits<XArray>::Type const& xMean
                            , typename hidden::Traits<YArray>::Type const& yMean
                            , bool unbiased = false
                          )
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(XArray);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(YArray);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
  typedef typename hidden::Traits<XArray>::Type XType;
  typedef typename hidden::Traits<YArray>::Type YType;
  typedef typename hidden::Promote<XType, YType>::result_type Type;

#ifdef STK_BOUNDS_CHECK
  if (X.range() != Y.range())
      STKRUNTIME_ERROR_NO_ARG(Error in Stat::covarianceWithFixedMeanSafe(X,Y,W,xMu,yMu,unbiased),ranges are not the sames);
#endif

  Type xsum  = 0.0, ysum = 0.0, xdev, ydev, sumWeights = 0.0, sum2Weights = 0.0, cov = 0.0;
  for (int i=X.begin(); i<X.end(); i++)
  {
    xsum += (xdev = X[i] - xMean);
    ysum += (ydev = Y[i] - yMean);
    Type Wi = std::abs((Real)W[i]);
    cov += Wi * (xdev*ydev);
    sumWeights  += Wi;
    sum2Weights += Wi * Wi;
  }
  // compute the variance
  return (unbiased) ?
       (sumWeights*sumWeights > sum2Weights) ?
         (cov - xsum*ysum)/(sumWeights - sum2Weights/sumWeights) : 0.
//         (cov - xsum*ysum/sumWeights)/(sumWeights - sum2Weights/sumWeights) : 0.
       : (sumWeights) ? (cov - xsum*ysum)/(sumWeights) : 0.;
}

/** @ingroup StatDesc
 *
 *  The covariance between X and Y will be computed using the formula
 *  @see Mark Galassi, Jim Davies, James Theiler,
 *  Brian Gough, Gerard Jungman, Michael Booth, and Fabrice Rossi.
 *  http://www.gnu.org/software/gsl/manual GNU Scientific Library -
 *   Reference manual, Version 1.15, 2011.
 *  http://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.html
 *  Sec. 21.7 Weighted Samples
 *
 *  \f$ cov(X,Y)=\frac{\sum_{i=1}^{N}w_i}{\left(\sum_{i=1}^{N}w_i\right)^2-\sum_{i=1}^{N}w_i^2}
 *  \sum_{i=1}^N w_i \left(  x_{i}-\mu_{x} \right)  \left( y_{i}-\mu_{y} \right). \f$
 *
 *  @param X,Y,W variables and weights
 *  @param xMu,yMu means of the variables
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template<class XArray, class YArray, class Weights>
Real covarianceWithFixedMeanSafe( ExprBase<XArray> const& X
                                , ExprBase<YArray> const& Y
                                , ExprBase<Weights> const& W
                                , typename hidden::Traits<XArray>::Type const& xMu
                                , typename hidden::Traits<YArray>::Type const& yMu
                                , bool unbiased = false
                                )
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(XArray);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(YArray);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
  typedef typename hidden::Traits<XArray>::Type XType;
  typedef typename hidden::Traits<YArray>::Type YType;
  typedef typename hidden::Promote<XType, YType>::result_type Type;

#ifdef STK_BOUNDS_CHECK
  if (X.range() != Y.range())
      STKRUNTIME_ERROR_NO_ARG(Error in Stat::covarianceWithFixedMeanSafe(X,Y,W,xMu,yMu,unbiased),ranges are not the sames);
#endif

  // compute covariance
  Type xsum  = 0.0, ysum = 0.0, xdev, ydev, sumWeights = 0.0, sum2Weights = 0.0, cov = 0.0;
  for (int i=X.begin(); i<X.end(); i++)
  {
    if  (Arithmetic<Type>::isFinite(X[i])
      && Arithmetic<Type>::isFinite(Y[i])
      && Arithmetic<Type>::isFinite(W[i])
      )
    {
      xsum += (xdev = X[i] - xMu); // deviation from the mean
      ysum += (ydev = Y[i] - yMu); // deviation from the mean
      Type Wi = std::abs((Real)W[i]);
      cov         += Wi * (xdev*ydev);         // cross product
      sumWeights  += Wi;           // sum absolute weights
      sum2Weights += Wi * Wi;     // sum squared weights
    }
  }
  // compute the variance
  return (unbiased) ?
       (sumWeights*sumWeights > sum2Weights) ?
         (cov - xsum*ysum)/(sumWeights - sum2Weights/sumWeights) : 0.
//         (cov - xsum*ysum/sumWeights)/(sumWeights - sum2Weights/sumWeights) : 0.
       : (sumWeights) ? (cov - xsum*ysum)/(sumWeights) : 0.;
}

/**  @ingroup StatDesc
 *  Compute the covariance matrix using the column of the data set V.
 *  \f[ \hat{\Sigma} = \frac{1}{n}
 *                       \sum_{i=1}^n (V_i-\hat{\mu}) (V_i-\hat{\mu})^T.
 *  \f]
 *  @param V variable
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template < class Array >
CArraySquare<typename Array::Type, Array::sizeCols_>
covariance( ExprBase<Array> const& V, bool unbiased = false)
{
  CArraySquare<typename Array::Type, Array::sizeCols_> cov_(V.cols());
  // typename Array::Row mean;
  typename hidden::FunctorTraits<Array, MeanOp>::Row mean;
  // compute the mean
  mean.move(Stat::mean(V.asDerived()));
  for (int j= cov_.begin(); j< cov_.end(); j++)
  {
    cov_(j, j) = varianceWithFixedMean(V.col(j), mean[j], unbiased);
    for (int i= cov_.begin(); i<j; i++)
    { cov_(j,i) = ( cov_(i, j) = covarianceWithFixedMean(V.col(i), V.col(j), mean[i], mean[j], unbiased));}
  }
  return cov_;
}

/**  @ingroup StatDesc
 *  Compute the covariance matrix using the rows of the data set V.
 *  \f[ \hat{\Sigma} = \frac{1}{p}
 *                       \sum_{j=1}^p (V_j-\hat{\mu})^T (V_j-\hat{\mu}).
 *  \f]
 *  @param V variable
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template < class Array >
CArraySquare<typename Array::Type, Array::sizeRows_>
covarianceByRow( ExprBase<Array> const& V, bool unbiased = false)
{
  CArraySquare<typename Array::Type, Array::sizeRows_> cov_(V.rows());
  typename hidden::FunctorTraits<Array, MeanOp>::Col mean;
  // compute the mean
  mean.move(Stat::meanByRow(V.asDerived()));
  for (int j= cov_.begin(); j< cov_.end(); j++)
  {
    cov_(j, j) = varianceWithFixedMean(V.row(j), mean[j], unbiased);
    for (int i= cov_.begin(); i<j; i++)
    { cov_(j,i) = ( cov_(i, j) = covarianceWithFixedMean(V.row(i), V.row(j), mean[i], mean[j], unbiased));}
  }
  return cov_;
}

/**  @ingroup StatDesc
 *  Compute the weighted covariance matrix using the column of the data set V.
 *  \f[ \hat{\sigma}^2 = \frac{1}{n}
 *                       \sum_{i=1}^n  w_{i} w_j (V_i-\hat{\mu}) (V_i-\hat{\mu})^T.
 *  \f]
 *  @param V,W the variables and the weights
 *  @param unbiased @c true if we want an unbiased estimate of the covariance,
 *  @c false otherwise (default is @c false)
 **/
template <class Array, class Weights >
CArraySquare<typename Array::Type, Array::sizeCols_>
covariance( ExprBase<Array> const& V, Weights const& W, bool unbiased = false)
{
  CArraySquare<typename Array::Type, Array::sizeCols_> cov_(V.cols());
  typename hidden::FunctorTraits<Array, MeanOp>::Row mean;
  // compute the mean
  mean.move(Stat::mean(V.asDerived(), W.asDerived()));
  for (int j= cov_.begin(); j< cov_.end(); j++)
  {
    cov_(j, j) = varianceWithFixedMean(V.col(j), W, mean[j], unbiased);
    for (int i= cov_.begin(); i<j; i++)
    { cov_(j,i) = ( cov_(i, j) = covarianceWithFixedMean(V.col(i), V.col(j), W, mean[i], mean[j], unbiased));}
  }
  return cov_;
}

/**  @ingroup StatDesc
 *  Compute the weighted covariance matrix using the rows of the data set V.
 *  \f[ \hat{\sigma}^2 = \frac{1}{n}
 *                       \sum_{i=1}^n  w_{i} w_j (V_i-\hat{\mu}) (V_i-\hat{\mu})^T.
 *  \f]
 *  @param V,W the variables and the weights
 *  @param unbiased @c true if we want an unbiased estimate of the covariance,
 *  @c false otherwise (default is @c false)
 **/
template <class Array, class Weights >
CArraySquare<typename Array::Type, Array::sizeRows_>
covarianceByRow( ExprBase<Array> const& V, Weights const& W, bool unbiased = false)
{
  CArraySquare<typename Array::Type, Array::sizeRows_> cov_(V.rows());
  typename hidden::FunctorTraits<Array, MeanOp>::Col mean;
  mean.move(Stat::meanByRow(V.asDerived(), W.asDerived()));
  for (int j= cov_.begin(); j< cov_.end(); j++)
  {
    cov_(j, j) = varianceWithFixedMean(V.row(j), W, mean[j], unbiased);
    for (int i= cov_.begin(); i<j; i++)
    { cov_(j,i) = ( cov_(i, j) = covarianceWithFixedMean(V.row(i), V.row(j), W, mean[i], mean[j], unbiased));}
  }
  return cov_;
}

/**  @ingroup StatDesc
 *  Compute the covariance matrix using the columns of the data set V.
 *  \f[ \hat{\Sigma} = \frac{1}{n}
 *                       \sum_{i=1}^n (V_i-mu) (V_i-mu)^T.
 *  \f]
 *  @param V variable
 *  @param mean mean of the variables
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template < class Array, class Mean >
CArraySquare<typename Array::Type, Array::sizeCols_>
covarianceWithFixedMean( ExprBase<Array> const& V, ExprBase<Mean> const& mean, bool unbiased = false)
{
#ifdef STK_BOUNDS_CHECK
  if (V.cols()!=mean.range()) STKRUNTIME_ERROR_NO_ARG(covarianceWithFixedMean,V.cols()!=mean.range());
#endif
  CArraySquare<typename Array::Type, Array::sizeCols_> cov_(V.cols());
  for (int j= cov_.begin(); j< cov_.end(); j++)
  {
    cov_(j, j) = varianceWithFixedMean(V.col(j), mean[j], unbiased);
    for (int i= cov_.begin(); i<j; i++)
    { cov_(j,i) = ( cov_(i, j) = covarianceWithFixedMean(V.col(i), V.col(j), mean[i], mean[j], unbiased));}
  }
  return cov_;
}

/**  @ingroup StatDesc
 *  Compute the covariance matrix using the rows of the data set V.
 *  \f[ \hat{\Sigma} = \frac{1}{n}
 *                       \sum_{i=1}^n (V_i-mu)^T (V_i-mu).
 *  \f]
 *  @param V variable
 *  @param mean mean of the variables
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template < class Array, class Mean >
CArraySquare<typename Array::Type, Array::sizeRows_>
covarianceWithFixedMeanByRow( ExprBase<Array> const& V, ExprBase<Mean> const& mean, bool unbiased = false)
{
#ifdef STK_BOUNDS_CHECK
  if (V.rows()!=mean.range()) STKRUNTIME_ERROR_NO_ARG(covarianceWithFixedMean,V.rows()!=mean.range());
#endif
  CArraySquare<typename Array::Type, Array::sizeRows_> cov_(V.rows());
  for (int j= cov_.begin(); j< cov_.end(); j++)
  {
    cov_(j, j) = varianceWithFixedMean(V.row(j), mean[j], unbiased);
    for (int i= cov_.begin(); i<j; i++)
    { cov_(j,i) = ( cov_(i, j) = covarianceWithFixedMean(V.row(i), V.row(j), mean[i], mean[j], unbiased));}
  }
  return cov_;
}

/**  @ingroup StatDesc
 *  Compute the weighted covariance matrix using the column of the data set V.
 *  \f[ \hat{\sigma}^2 = \frac{1}{n}
 *                       \sum_{i=1}^n  w_{i} w_j (V_i-\hat{\mu}) (V_i-\hat{\mu})^T.
 *  \f]
 *  @param V,W the variables and the weights
 *  @param mean the (weighted) mean of the variables
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template <class Array, class Weights, class Mean >
CArraySquare<typename Array::Type, Array::sizeCols_>
covarianceWithFixedMean(  ExprBase<Array> const& V, Weights const& W, ExprBase<Mean> const& mean, bool unbiased = false)
{
#ifdef STK_BOUNDS_CHECK
  if (V.cols()!=mean.range()) STKRUNTIME_ERROR_NO_ARG(covarianceWithFixedMean,V.cols()!=mean.range());
  if (W.range()!=mean.range()) STKRUNTIME_ERROR_NO_ARG(covarianceWithFixedMean,W.range()!=mean.range());
#endif

  CArraySquare<typename Array::Type, Array::sizeCols_> cov_(V.cols());
  for (int j= cov_.begin(); j< cov_.end(); j++)
  {
    cov_(j, j) = varianceWithFixedMean(V.col(j), W, mean[j], unbiased);
    for (int i= cov_.begin(); i<j; i++)
    { cov_(j,i) = ( cov_(i, j) = covarianceWithFixedMean(V.col(i), V.col(j), W, mean[i], mean[j], unbiased));}
  }
  return cov_;
}

/**  @ingroup StatDesc
 *  Compute the weighted covariance matrix using the rows of the data set V.
 *  \f[ \hat{\sigma}^2 = \frac{1}{n}
 *                       \sum_{i=1}^n  w_{i} w_j (V_i-\hat{\mu})^T (V_i-\hat{\mu}).
 *  \f]
 *  @param V,W the variables and the weights
 *  @param mean the (weighted) mean of the variables
 *  @param unbiased @c true if we want an unbiased estimate of the variance,
 *  @c false otherwise (default is @c false)
 **/
template <class Array, class Weights, class Mean >
CArraySquare<typename Array::Type, Array::sizeRows_>
covarianceWithFixedMeanByRow(  ExprBase<Array> const& V, Weights const& W, ExprBase<Mean> const& mean, bool unbiased = false)
{
#ifdef STK_BOUNDS_CHECK
  if (V.rows()!=mean.range()) STKRUNTIME_ERROR_NO_ARG(covarianceWithFixedMean,V.cols()!=mean.range());
  if (W.range()!=mean.range()) STKRUNTIME_ERROR_NO_ARG(covarianceWithFixedMean,W.range()!=mean.range());
#endif

  CArraySquare<typename Array::Type, Array::sizeRows_> cov_(V.rows());
  for (int j= cov_.begin(); j< cov_.end(); j++)
  {
    cov_(j, j) = varianceWithFixedMean(V.row(j), W, mean[j], unbiased);
    for (int i= cov_.begin(); i<j; i++)
    { cov_(j,i) = ( cov_(i, j) = covarianceWithFixedMean(V.row(i), V.row(j), W, mean[i], mean[j], unbiased));}
  }
  return cov_;
}

}  // namespace Stat

}  // namespace STK

#endif /*STK_STAT_COVARIANCE_H */
