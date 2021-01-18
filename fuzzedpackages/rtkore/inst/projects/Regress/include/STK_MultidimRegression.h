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
 * created on: 27 oct. 2010
 * Purpose:  .
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_MultidimRegression.h
 *  @brief In this file we define the MultidimRegression class.
 **/

#ifndef STK_MULTIDIMREGRESSION_H
#define STK_MULTIDIMREGRESSION_H

#include <Arrays/include/STK_Array2D.h> // for coefs
#include <Algebra/include/STK_InvertMatrix.h>
#include "STK_IRegression.h"

namespace STK
{

/** @brief The @c MultidimRegression class allows to regress a multidimensional
 *  output variable among a multivariate explanation variable.
 */
template<class Array, class Weight>
class MultidimRegression: public IRegression<Array, Array, Weight>
{
  public:
    typedef IRegression<Array, Array, Weight> Base;
    using Base::p_x_;
    using Base::p_y_;
    /** Constructor.
     * @param y,x Variates to predict and co-variates
     */
    MultidimRegression( Array const* y =0, Array const* x =0);
    /** Destructor. */
    virtual ~MultidimRegression() {}
    /** @return the coefficients */
    inline Array const& coefs() const { return coefs_;}
    /** @return the extrapolated values y from the value @c x.
     *  Given the data set @c x will compute the values \f$ y = x.\hat{\beta} \f$.
     *  The coefficients @c coefs_ have to be estimated previously.
     *  @param x the input data set
     */
    virtual Array extrapolate(Array const& x) const;

  protected:
    ArrayXX coefs_;

  private:
    /** compute the regression function. */
    virtual bool regressionStep();
    /** compute the weighted regression function.
     * @param weights the weights of the samples
     **/
    virtual bool regressionStep(Weight const& weights);
    /** Compute the predicted outputs by the regression function. */
    virtual bool predictionStep();
    /** Compute the number of parameter of the regression function.
     * @return the number of parameter of the regression function
     **/
    inline virtual int computeNbFreeParameter() const
    { return coefs_.sizeCols() * coefs_.sizeRows(); }
};

template<class Array, class Weight>
MultidimRegression<Array,Weight>::MultidimRegression( Array const* y, Array const* x)
                                                   : Base(y, x)
                                                    , coefs_()
{}

/* compute the regression function. */
template<class Array, class Weight>
bool MultidimRegression<Array,Weight>::regressionStep()
{
  // compute X'X
  ArraySquareX prod;
  prod.move(multLeftTranspose(p_x_->asDerived()));
  // compute (X'X)^{-1}
//  GInvertSymMatrix<ArraySquareX>()(prod);
  // compute (X'X)^{-1}X'Y
  coefs_.move(mult(invert(prod.symmetrize()), multLeftTranspose(p_x_->asDerived(), p_y_->asDerived())));
  //coefs_.move(mult(prod, multLeftTranspose(p_x_->asDerived(), p_y_->asDerived())));
  return true;
}

/* compute the regression function. */
template<class Array, class Weight>
bool MultidimRegression<Array,Weight>::regressionStep(Weight const& weights)
{
  // compute X'WX
  ArraySquareX prod;
  prod.move(weightedMultLeftTranspose(p_x_->asDerived(), weights));
  // compute (X'WX)^{-1}
  //GInvertSymMatrix<ArraySquareX>()(prod);
  // compute (X'WX)^{-1}X'WY
  coefs_.move(mult(invert(prod.symmetrize()), wmultLeftTranspose(p_x_->asDerived(), p_y_->asDerived(), weights)));
//  coefs_.move(mult(prod, wmultLeftTranspose(p_x_->asDerived(), p_y_->asDerived(), weights)));
  return true;
}

/* Compute the predicted outputs by the regression function. */
template<class Array, class Weight>
bool MultidimRegression<Array,Weight>::predictionStep()
{ this->predicted_.move(mult(*p_x_, coefs_));  return true;}

/* @brief Extrapolate the the values @c y from the value @c x.
 *  Given the data set @c x will compute the values \f$ y = \hat{f}(x) \f$.
 *  The regression function @e f have to be estimated previously.
 */
template<class Array, class Weight>
Array MultidimRegression<Array,Weight>::extrapolate( Array const& x) const
{ return(mult(x, coefs_));}

} // namespace STK

#endif /* STK_MULTIDIMREGRESSION_H */
