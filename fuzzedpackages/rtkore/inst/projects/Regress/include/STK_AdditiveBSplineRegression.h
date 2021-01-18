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
 * Purpose: definition of the AdditiveBSplineRegression class.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_AdditiveBSplineRegression.h
 *  @brief In this file we define the AdditiveBSplineRegression class.
 **/

#ifndef STK_ADDITIVEBSPLINEREGRESSION_H
#define STK_ADDITIVEBSPLINEREGRESSION_H

#include "STK_AdditiveBSplineCoefficients.h"
#include "STK_IRegression.h"

#ifdef STKUSELAPACK
#include <Algebra/include/STK_lapack_MultiLeastSquare.h>
#else
#include <Algebra/include/STK_MultiLeastSquare.h>
#endif

namespace STK
{

/** @ingroup Regress
 *  @brief Compute an additive BSpline, multivalued, regression function using
 *  BSpline basis.
 */
template <class YArray, class XArray, class Weights = VectorX>
class AdditiveBSplineRegression: public IRegression<YArray, XArray, Weights>
{
  public:
    typedef IRegression<YArray, XArray, Weights> Base;
    typedef Regress::KnotsPosition KnotsPosition;
    using Base::p_x_;
    using Base::p_y_;
    using Base::predicted_;
    using Base::residuals_;
    /** Constructor.
     * @param p_y p-dimensional array of output to fit
     * @param p_x d-dimensional array of predictor
     * @param nbControlPoints number of control points of the spline
     * @param degree degree of the BSpline basis
     * @param position position of the knots to used
     **/
    AdditiveBSplineRegression( YArray const* p_y
                             , XArray const* p_x
                             , int nbControlPoints
                             , int degree = 3
                             , KnotsPosition const& position = Regress::uniformKnotsPositions_
                             );
    /** Constructor.
     * @param y p-dimensional array of output to fit
     * @param x d-dimensional array of predictor
     * @param nbControlPoints number of control points of the spline
     * @param degree degree of the BSpline basis
     * @param position position of the knots to used
     **/
    AdditiveBSplineRegression( YArray const& y
                             , XArray const& x
                             , int nbControlPoints
                             , int degree = 3
                             , KnotsPosition const& position = Regress::uniformKnotsPositions_
                             );
    /** virtual destructor. */
    virtual ~AdditiveBSplineRegression() {}
    /** @return the degree of the B-spline curves */
    inline int degree() const { return degree_;}
    /** @return the number of control points of the B-spline curves */
    inline int nbControlPoints() const { return nbControlPoints_;}
    /** @return the control points of the B-spline curves */
    inline YArray const& controlPoints() const { return controlPoints_; }
    /** This is a matrix of size (p_x_->range(), 0:lastControlPoints).
     *  @return the coefficients of the B-spline curves
     **/
    inline ArrayXX const& coefficients() const { return coefs_.coefficients();}
    /** @return The extrapolated values of y from the value @c x.
     *  Given the data set @c x will compute the values \f$ y = \psi(x) \hat{\beta} \f$
     *  where \f$ \psi \f$ represents the B-spline basis functions and \f$ \hat{beta} \f$
     *  the estimated coefficients.
     *  @param x the input data set
     **/
    virtual YArray extrapolate( XArray const& x) const;

  protected:
    /** number of control points of the B-spline curve. */
    int nbControlPoints_;
    /** degree of the B_Spline curve */
    int degree_;
    /** method of position of the knots of the B-spline curve */
    KnotsPosition position_;
    /** Coefficients of the regression matrix */
    AdditiveBSplineCoefficients<XArray> coefs_;
    /** Estimated control points of the B-spline curve */
    YArray controlPoints_;
    /** compute the coefficients of the BSpline basis. This method will be
     * called in the base class @c IRegression::run()
     **/
    virtual bool initializeStep();
    /** Compute the regression function. */
    virtual bool regressionStep();
    /** Compute the weighted regression.
     *  @param weights weights of the samples
     **/
    virtual bool regressionStep(Weights const& weights);
    /** Compute the predicted outputs. */
    virtual bool predictionStep();
    /** @return the number of parameter of the regression function */
    inline virtual int computeNbFreeParameter() const
    { return controlPoints_.sizeCols() * controlPoints_.sizeRows(); }
};

/* Constructor.
 * @param p_y p-dimensional array of output to fit
 * @param p_x d-dimensional array of predictor
 * @param nbControlPoints number of control points of the spline
 * @param degree degree of the BSpline basis
 * @param position position of the knots to used
 **/
template <class YArray, class XArray, class Weights>
AdditiveBSplineRegression< YArray, XArray, Weights>::AdditiveBSplineRegression(
                                                      YArray const* p_y
                                                    , XArray const* p_x
                                                    , int nbControlPoints
                                                    , int degree
                                                    , KnotsPosition const& position
                                                    )
                                                   : Base(p_y, p_x)
                                                    , nbControlPoints_(nbControlPoints)
                                                    , degree_(degree)
                                                    , position_(position)
                                                    , coefs_(p_x, nbControlPoints_, degree_, position_)
                                                    , controlPoints_()
{}
/* Constructor.
 * @param y p-dimensional array of output to fit
 * @param x d-dimensional array of predictor
 * @param nbControlPoints number of control points of the spline
 * @param degree degree of the BSpline basis
 * @param position position of the knots to used
 **/
template <class YArray, class XArray, class Weights>
AdditiveBSplineRegression< YArray, XArray, Weights>::AdditiveBSplineRegression( YArray const& y
                                                    , XArray const& x
                                                    , int nbControlPoints
                                                    , int degree
                                                    , KnotsPosition const& position
                                                    )
                                                   : Base(&y, &x)
                                                    , nbControlPoints_(nbControlPoints)
                                                    , degree_(degree)
                                                    , position_(position)
                                                    , coefs_(&x, nbControlPoints_, degree_, position_)
                                                    , controlPoints_()
{}

/* compute the coefficients of the BSpline basis. This method will be
 * called in the base class @c IRegression::run()
 **/
template <class YArray, class XArray, class Weights>
bool AdditiveBSplineRegression< YArray, XArray, Weights>::initializeStep()
{
  coefs_.setData(p_x_, nbControlPoints_, degree_, position_);
  if (!coefs_.run())
  {
    this->msg_error_ = coefs_.error();
    return false;
  }
  return true;
}

/* compute the regression function. */
template <class YArray, class XArray, class Weights>
bool AdditiveBSplineRegression< YArray, XArray, Weights>::regressionStep()
{
  // coefs_.coefficients() is Array2D
#ifdef STKUSELAPACK
  lapack::MultiLeastSquare<YArray, ArrayXX> reg(*p_y_, coefs_.coefficients(), false, false);
#else
  MultiLeastSquare<YArray, ArrayXX> reg(*p_y_, coefs_.coefficients(), false, false);
#endif
  if (!reg.run())
  {
    this->msg_error_ = reg.error();
    return false;
  }
  controlPoints_.move(reg.x());
  return true;
}

/* compute the regression function. */
template <class YArray, class XArray, class Weights>
bool AdditiveBSplineRegression< YArray, XArray, Weights>::regressionStep(Weights const& weights)
{
#ifdef STKUSELAPACK
  lapack::MultiLeastSquare<YArray, ArrayXX> reg(*p_y_, coefs_.coefficients(), false, false);
#else
  MultiLeastSquare<YArray, ArrayXX> reg(*p_y_, coefs_.coefficients(), false, false);
#endif
  if (!reg.run(weights))
  {
    this->msg_error_ = reg.error();
    return false;
  }
  controlPoints_ = reg.x();
  return true;
}

/* Compute the predicted outputs by the regression function. */
template <class YArray, class XArray, class Weights>
bool AdditiveBSplineRegression< YArray, XArray, Weights>::predictionStep()
{
  predicted_ = coefs_.coefficients() * controlPoints_;
  return true;
}


/* @brief Extrapolate the values @c y from the value @c x.
 *  Given the data set @c x will compute the values \f$ y = \psi(x) \hat{\beta} \f$
 *  where \f$ \psi \f$ represents the B-spline basis functions and \f$ \hat{beta} \f$
 *  the estimated coefficients.
 */
template <class YArray, class XArray, class Weights>
YArray AdditiveBSplineRegression< YArray, XArray, Weights>::extrapolate( XArray const& x) const
{
  YArray res = coefs_.extrapolate(x) * controlPoints_;
  return res;
}

} // namespace STK

#endif /* STK_ADDITIVEBSPLINEREGRESSION_H */
