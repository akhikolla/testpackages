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
 * created on: 31 juil. 2010
 * Purpose: definition of the BsplineRegression class.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_BSplineRegression.h
 *  @brief In this file we define the BsplineRegression class.
 **/

#ifndef STK_BSPLINEREGRESSION_H
#define STK_BSPLINEREGRESSION_H

#include "STK_Regress_Util.h"
#include "STK_BSplineCoefficients.h"
#include "STK_IRegression.h"

#include <Algebra/include/STK_InvertMatrix.h>

namespace STK
{

/** @brief Compute a BSpline, multi-valued, regression function using BSpline
 *  basis.
 */
template <class YArray, class XVector, class Weights = VectorX>
class BSplineRegression: public IRegression<YArray, XVector, Weights>
{
  public:
    typedef IRegression<YArray, XVector, Weights> Base;
    typedef Regress::KnotsPosition KnotsPosition;
    using Base::p_x_;
    using Base::p_y_;
    using Base::predicted_;
    using Base::residuals_;
    /** Constructor.
     * @param p_y d-dimensional array of output to fit
     * @param p_x uni-dimensional array of predictor
     * @param nbControlPoints number of control points of the spline
     * @param degree degree of the BSpline basis
     * @param position position of the knots to used
     **/
    BSplineRegression( YArray const* p_y
                     , XVector const* p_x
                     , int const& nbControlPoints
                     , int const& degree = 3
                     , KnotsPosition const& position = Regress::uniformKnotsPositions_
                     );
    /** virtual destructor. */
    inline virtual ~BSplineRegression();
    /** @return the degree of the B-spline curve */
    inline int degree() const { return degree_;}
    /** @return the number of control points of the B-spline curve */
    inline int nbControlPoints() const { return nbControlPoints_;}
    /** @return the control points of the B-spline curve */
    inline YArray const& controlPoints() const { return controlPoints_; }
    /**  @return the knots of the B-spline curve */
    inline VectorX const& knots() const { return coefs_.knots(); }
    /** @return the coefficients of the B-spline curve */
    inline YArray const& coefficients() const { return coefs_.coefficients();}
    /** @return the Extrapolates values of y from the value @c x.
     *  Given the data set @c x will compute the values \f$ y = \psi(x) \hat{\beta} \f$
     *  where \f$ \psi \f$ represents the B-spline basis functions and \f$ \hat{beta} \f$
     *  the estimated coefficients.
     *  @param x the input data set
     */
    virtual YArray extrapolate( XVector const& x) const;

  protected:
    /** number of control points of the B-spline curve. */
    int nbControlPoints_;
    /** degree of the B_Spline curve */
    int degree_;
    /** method of position of the knots of the B-spline curve */
    KnotsPosition position_;
    /** Coefficients of the regression matrix */
    BSplineCoefficients<XVector> coefs_;
    /** Estimated control points of the B-spline curve */
    YArray controlPoints_;
    /** Compute the coefficients of the BSpline basis. This method is triggered
     *  by the base class @c IRegression::run()
     **/
    virtual bool initializeStep();
    /** Compute the regression function. This method is triggered by
     *  the base class @c IRegression::run() after initializeStep()
     **/
    virtual bool regressionStep();
    /** Compute the regression function. This method is triggered
     *  by the base class @c IRegression::run(weights) after initializeStep()
     *  @param weights the weights of the samples
     **/
    virtual bool regressionStep(Weights const& weights);
    /** Compute the predicted outputs by the regression function. This method
     *  is triggered by the base class @c IRegression::run() after initializeStep()
     **/
    virtual bool predictionStep();
    /** Compute the number of parameter of the regression function.
     * @return the number of parameter of the regression function
     **/
    inline virtual int computeNbFreeParameter() const
    { return controlPoints_.sizeCols() * controlPoints_.sizeRows(); }
};

template <class YArray, class XVector, class Weights>
BSplineRegression<YArray, XVector, Weights>::BSplineRegression( YArray const* p_y
                                                              , XVector const* p_x
                                                              , int const& nbControlPoints
                                                              , int const& degree
                                                              , const KnotsPosition& position
                                                              )
                                                             : Base(p_y, p_x)
                                                              , nbControlPoints_(nbControlPoints)
                                                              , degree_(degree)
                                                              , position_(position)
                                                              , coefs_(*p_x, nbControlPoints_, degree_, position_)
                                                              , controlPoints_()
{ }

template <class YArray, class XVector, class Weights>
BSplineRegression<YArray, XVector, Weights>::~BSplineRegression()
{}

template <class YArray, class XVector, class Weights>
bool BSplineRegression<YArray, XVector, Weights>::initializeStep()
{ return coefs_.run();}
/* compute the regression function. */
template <class YArray, class XVector, class Weights>
bool BSplineRegression<YArray, XVector, Weights>::regressionStep()
{
  // compute X'X
  ArraySquareX prod = coefs_.coefficients().transpose() * coefs_.coefficients();
  // compute (X'X)^{-1}
//  GInvertSymMatrix<ArraySquareX> inv;
//  inv(prod);

  // compute (X'X)^{-1}X'Y
  controlPoints_ = invert(prod.symmetrize()) * (coefs_.coefficients().transpose() * p_y_->asDerived());
  return true;
}

/* compute the regression function. */
template <class YArray, class XVector, class Weights>
bool BSplineRegression<YArray, XVector, Weights>::regressionStep(Weights const& weights)
{
  // compute X'X
  ArraySquareX prod = coefs_.coefficients().transpose() * weights.diagonalize() * coefs_.coefficients();
  // compute (X'X)^{-1}
  // GInvertSymMatrix<ArraySquareX> inv;
  // inv(prod);
  // compute (X'X)^{-1}X'Y
  controlPoints_ = invert(prod.symmetrize()) * coefs_.coefficients().transpose() * weights.diagonalize() * p_y_->asDerived();
  //controlPoints_ = prod * coefs_.coefficients().transpose() * weights.diagonalize() * p_y_->asDerived();
  return true;
}

/* Compute the predicted outputs by the regression function. */
template <class YArray, class XVector, class Weights>
bool BSplineRegression<YArray, XVector, Weights>::predictionStep()
{
  predicted_ = coefs_.coefficients() * controlPoints_;
  return true;
}

/* @brief Extrapolate the values @c y from the value @c x.
 *  Given the data set @c x will compute the values \f$ y = x.\hat{\beta} \f$.
 *  The coefficients @c coefs_ have to be estimated previously.
 *  @param x the input data set
 *  @param y the output (extrapolated) data set
 */
template <class YArray, class XVector, class Weights>
YArray BSplineRegression<YArray, XVector, Weights>::extrapolate( XVector const& x) const
{
  YArray res = coefs_.extrapolate(x) * controlPoints_;
  return res;
}

} // namespace STK

#endif /* STK_BSPLINEREGRESSION_H */
