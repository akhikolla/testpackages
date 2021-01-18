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
 * created on: 23 juin 2010
 * Purpose:  Interface base class for regression methods.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IRegression.h
 *  @brief In this file we define the Interface base class IRegression.
 **/

#ifndef STK_IREGRESSION_H
#define STK_IREGRESSION_H

#include <STKernel/include/STK_String.h>

namespace STK
{

/** @ingroup Regress
 * @brief Interface base class for Regression methods.
 *
 * Regression models involve the following variables:
 * - The <em>unknown parameters</em>, denoted as <em> β</em>, which may represent
 *  a scalar or a vector.
 * - The <em>independent variables</em>, <em>X</em>.
 * - The <em>dependent variable</em>, <em>Y</em>.
 *
 * In various fields of application, different terminologies are used in place
 * of dependent and independent variables. A regression model relates <em>Y</em>
 * to a function of <em>X</em> and <em>β</em>.
 * \f$ Y \approx f (\mathbf {X}, \boldsymbol{\beta} ) \f$
 *
 * The approximation is usually formalized as \f$ E(Y|X)=f(X,β) \f$. To carry out
 * regression analysis, the form of the function <em>f</em> must be specified.
 * Sometimes the form of this function is based on knowledge about the relationship
 * between <em>Y</em> and <em>X</em> that does not rely on the data. If no such
 *  knowledge is available, a flexible or convenient form for <em>f</em> is chosen.
 *
 * Assume now that the vector of unknown parameters <em>β</em> is of length <em>k</em>.
 * In order to perform a regression analysis the user must provide information
 * about the dependent variable <em>Y</em>:
 * - If <em>N</em> data points of the form (<em>Y</em>,<em>X</em>) are observed,
 * where <em>N</em> < <em>k</em>, most classical approaches to regression analysis
 * cannot be performed: since the system of equations defining the regression model
 * is under-determined, there is not enough data to recover <em>β</em>.
 * If exactly <em>N</em>=<em>k</em> data points are observed, and the function
 * <em>f</em> is linear, the equations <em>Y</em>=<em>f</em>(<em>X</em>, <em>β</em>)
 * can be solved exactly rather than approximately. This reduces to solving a set
 * of <em>N</em> equations with <em>N</em> unknowns (the elements of <em>β</em>),
 * which has a unique solution as long as the <em>X</em> are linearly independent.
 * If <em>f</em> is nonlinear, a solution may not exist, or many solutions may exist.
 * - The most common situation is where <em>N</em> > <em>k</em> data points are
 * observed. In this case, there is enough information in the data to estimate a
 * unique value for <em>β</em> that best fits the data in some sense, and the
 * regression model when applied to the data can be viewed as an overdetermined
 * system in <em>β</em>.
 *
 * In the last case, the regression analysis provides the tools for:
 * #- Finding a solution for unknown parameters <em>β</em> that will, for example,
 * minimize the distance between the measured and predicted values of the
 * dependent variable <em>Y</em> (also known as method of least squares).
 * #- Under certain statistical assumptions, the regression analysis uses the
 * surplus of information to provide statistical information about the
 * unknown parameters <em>β</em> and predicted values of the dependent variable
 * <em>Y</em>.
 *
 * In this interface, the pure virtual function to implement are
 * @code
 *   virtual bool regressionStep()=0;
 *   virtual bool regressionStep( Weights const& weights) =0;
 *   virtual bool predictionStep() =0;
 *   virtual int computeNbFreeParameter() const =0;
 *   virtual bool extrapolate( XArray const& x, YArray& y) const =0;
 * @endcode
 * The virtual function
 * @code
 *   virtual void initializeStep();
 * @endcode
 * can be overloaded.
 */
template <class YArray, class XArray, class Weights>
class IRegression: public IRunnerSupervised<YArray, XArray, Weights>
{
  protected:
    typedef IRunnerSupervised<YArray, XArray, Weights> Base;
    using Base::p_x_;
    using Base::p_y_;
    /** Default constructor. Initialize the data members. */
    IRegression(): Base(), predicted_(), residuals_(), nbFreeParameter_(0) {}
    /** constructor
     *  @param p_y,p_x pointers on the y and x data sets
     **/
    IRegression( ArrayBase<YArray> const* p_y, ArrayBase<XArray> const* p_x)
              : Base((p_y == 0) ? 0 : p_y->asPtrDerived(), (p_x == 0) ? 0 : p_x->asPtrDerived())
               , predicted_()
               , residuals_()
               , nbFreeParameter_(0)
     {}
    /** Constructor. Initialize the data members.
     *  @param y,x arrays with the observed output and inputs of the model
     */
    IRegression( ArrayBase<YArray> const& y, ArrayBase<XArray> const& x)
              : Base(y.asDerived(), x.asDerived())
               , predicted_()
               , residuals_()
               , nbFreeParameter_(0)
     {}

  public:
    /** virtual destructor. */
    virtual ~IRegression() {}

    /** @return the predicted values */
    inline YArray const& predicted() const { return predicted_;}
    /** @return the residuals */
    inline YArray const& residuals() const {  return residuals_;}

    /** @return the pointer on the predicted values */
    inline YArray* p_predicted() { return &predicted_;}
    /**  @return the pointer on the residuals */
    inline YArray* p_residuals() { return &residuals_;}

    /**  @return the number of parameter of the regression function */
    inline int nbFreeParameter() const { return nbFreeParameter_;}

    /** run the computations. Default Implementation. */
    virtual bool run();
    /** run the weighted computations.
     *  @param weights weights of the samples
     **/
    virtual bool run( Weights const& weights);
    /** @return the compute values of @c y from the value @c x using the model.
     *  Given the data set @c x will compute the values \f$ y = \hat{f}(x) \f$.
     *  The regression function @e f has to be estimated previously.
     *  @param x the input data set
     */
    virtual YArray extrapolate(XArray const& x) const =0;

  protected:
    /** @brief perform any computation needed before the call of the regression
     *  method. Default implementation is do nothing.
     */
    virtual bool initializeStep() { return true;}
    /** @brief perform any computation needed after the call of the regression
     *  method. Default implementation is do nothing.
     */
    virtual bool finalizeStep() { return true;}
    /** @brief Compute the residuals of the model.
     * The residuals of the model are computed by computing the difference
     * between the observed outputs and the predicted outputs of the model.
     */
    inline bool residualsStep()
    {
      residuals_ = *p_y_ - predicted_;
      return true;
    }

    /** Container of the predicted output. */
    YArray predicted_;
    /** Container of the residuals. */
    YArray residuals_;

  private:
    /** number of parameter of the regression method. */
    int nbFreeParameter_;
    /** compute the regression function. */
    virtual bool regressionStep() =0;
    /** compute the weighted regression function.
     * @param weights the weights of the samples
     **/
    virtual bool regressionStep(Weights const& weights) =0;
    /** Compute the predicted outputs by the regression function and store the
     * result in the p_predicted_ array. */
    virtual bool predictionStep() =0;
    /** Compute the number of parameter of the regression function.
     * @return the number of parameter of the regression function
     **/
    virtual int computeNbFreeParameter() const =0;
};

/** run the computations. Default Implementation. */
template <class YArray, class XArray, class Weights>
bool IRegression<YArray,XArray,Weights>::run()
{
  // perform any initialization step needed before the regression step
  if (!initializeStep()) { return false;}
  // compute the regression
  if (!regressionStep()) { return false;}
  // Compute the number of parameter of the regression function.
  nbFreeParameter_ = computeNbFreeParameter();
  // compute predictions
  predictionStep();
  // compute residuals
  residualsStep();
  // perform any post-operation needed before the regression step
  finalizeStep();
  // return the result of the computations
  this->hasRun_ = true;
  return true;
}

/** run the computations. Default Implementation. */
template <class YArray, class XArray, class Weights>
bool IRegression<YArray,XArray,Weights>::run( Weights const& weights)
{
  // perform any pre-operation needed before the regression step
  if (!initializeStep()) { return false;}
  // compute weighted regression
  if (!regressionStep(weights)) { return false;}
  // Compute the number of parameter of the regression function.
  nbFreeParameter_ = computeNbFreeParameter();
  // create array of the predicted value and compute prediction
  predictionStep();
  // create array of the residuals and compute them
  residualsStep();
  // perform any post-operation needed before the regression step
  finalizeStep();
  // return the result of the computations
  this->hasRun_ = true;
  return true;
}

} // namespace STK

#endif /* STK_IREGRESSION_H */
