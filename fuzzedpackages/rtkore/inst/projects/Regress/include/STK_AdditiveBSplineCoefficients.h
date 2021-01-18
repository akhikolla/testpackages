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
 * created on: 25 juin 2010
 * Purpose:  Compute the coefficients of an additive B-spline manifold.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_AdditiveBSplineCoefficients.h
 *  @brief In this file we define the AdditiveBSplineCoefficients class.
 **/

#ifndef STK_ADDITIVEBSPLINECOEFFICIENTS_H
#define STK_ADDITIVEBSPLINECOEFFICIENTS_H

#include "STK_BSplineCoefficients.h"

namespace STK
{

/** @ingroup Regress
 *  @brief Compute the regression splines coefficients of an additive model.
 *
 * The method is described in @ref BSplineCoefficients documentation class and
 * repeated for each variables of the model. The number of control points, the
 * degree and the position of the knots are the same for all variables.
 *
 * If the input data set is a matrix of size (n,p) then the output matrix of
 * the coefficients @c Coefficients() is a matrix of size
 * (n, p*nbControlPoints) where p is the number of variables.
 */
template<class Array>
class AdditiveBSplineCoefficients: public IRunnerBase
{
  public:
    typedef typename hidden::Traits<Array>::Row RowVector;
    typedef typename hidden::Traits<Array>::Col ColVector;
    /** Constructor : initialize the data members. The number of knots is given
     *  by the formula nbKnots = nbControlPoints + degree +1.
     *  @param p_data a pointer on the input data values
     *  @param nbControlPoints number of control points
     *  @param degree degree of the B-spline curves
     *  @param position method to use for positioning the knots
     **/
    AdditiveBSplineCoefficients( Array const* p_data
                               , int nbControlPoints
                               , int degree = 3
                               , Regress::KnotsPosition position = Regress::uniformKnotsPositions_
                               )
                               : IRunnerBase()
                               , p_data_(p_data)
                               , nbKnots_(nbControlPoints + degree +1)
                               , nbControlPoints_(nbControlPoints)
                               , degree_(degree)
                               , position_(position)
                               , coefficients_()
    {}
    /** Constructor : initialize the data members. The number of knots is given
     *  by the formula nbKnots = nbControlPoints + degree +1.
     *  @param data the input data values
     *  @param nbControlPoints number of control points
     *  @param degree degree of the B-spline curves
     *  @param position method to use for positioning the knots
     **/
    AdditiveBSplineCoefficients( Array const& data
                               , int nbControlPoints
                               , int degree = 3
                               , Regress::KnotsPosition position = Regress::uniformKnotsPositions_
                               )
                               : IRunnerBase()
                               , p_data_(&data)
                               , nbKnots_(nbControlPoints + degree +1)
                               , nbControlPoints_(nbControlPoints)
                               , degree_(degree)
                               , position_(position)
                               , coefficients_()
    {}
    /** Copy constructor.
     *  @param coefs the coefficients to copy
     **/
    AdditiveBSplineCoefficients( AdditiveBSplineCoefficients const& coefs)
                               : IRunnerBase(coefs)
                               , p_data_(coefs.p_data_)
                               , nbKnots_(coefs.nbKnots_)
                               , nbControlPoints_(coefs.nbControlPoints_)
                               , degree_(coefs.degree_)
                               , position_(coefs.position_)
                               , coefficients_(coefs.coefficients_)
    {}
    /** Destructor. */
    virtual ~AdditiveBSplineCoefficients() {}
    /** clone pattern implementation */
    AdditiveBSplineCoefficients* clone() const { return new AdditiveBSplineCoefficients(*this);}
    /** Compute the coefficients of the B-spline curve for the given values.
     *  @param p_data the input data values
     *  @param nbControlPoints number of control points
     *  @param degree degree of the B-spline curves
     *  @param position method to use for positioning the knots
     **/
    void setData( Array const* p_data
                , int nbControlPoints
                , int degree = 3
                , Regress::KnotsPosition position = Regress::uniformKnotsPositions_
                );

    /** run the computations. */
    bool run();

    /** give the degree of the B-spline curves. */
    inline int degree() const { return degree_;}
    /** give the number of knots of the B-spline curves. */
    inline int nbKnots() const { return nbKnots_;}
    /** give the number of control points of the B-spline curves. */
    inline int nbControlPoints() const { return nbControlPoints_;}
    /** give the computed coefficients of the B-spline curves.
     *  This is a matrix of size (p_data_->range(), 0:lastControlPoints).
     **/
    inline ArrayXX const& coefficients() const { return coefficients_;}
    /** @return the extrapolated coefficients for the x matrix.
     *  @param x the values to extrapolate
     **/
    ArrayXX extrapolate(Array const& x) const;

  protected:
    /** the input data set */
    Array const* p_data_;
    /** number of knots of the B-spline curves.*/
    int nbKnots_;
    /** number of control points of the B-spline curves.*/
    int nbControlPoints_;
    /** degree of the B-splines curves. */
    int degree_;
    /** Method used in order to position the knots. */
    Regress::KnotsPosition position_;
    /** Array with the knots and coefficients in each dimension. */
    Array1D<BSplineCoefficients<ColVector> > coefs_;
    /** Array of the coefficients */
    ArrayXX coefficients_;
};

/* run the computations. */
template<class Array>
bool AdditiveBSplineCoefficients<Array>::run()
{
#ifdef STK_REGRESS_VERBOSE
  stk_cout << _T("In AdditiveBSplineCoefficients::run()\n");
#endif
  // check if there exists data
  if (!p_data_)
  {
    msg_error_ = STKERROR_NO_ARG(AdditiveBSplineCoefficients::run,data is not set);
    return false;
  }
  try
  {
    // resize the Array with the marginal coefficients
    coefs_.clear(); // necesssary as coefficient_ is a reference at second call
    coefs_.resize(p_data_->cols());

    // resize the matrix of all coefficient
    coefficients_.resize(p_data_->rows(), Range());
    // get dimensions
    for (int i=p_data_->beginCols(); i<p_data_->endCols(); i++)
    {
      // create a reference on the ith column of the data
      ColVector colData(p_data_->col(i), true);
      // set data to the i-th coefficient
      // WARNING: colData will be invalidate
      coefs_[i].setData(colData, nbControlPoints_, degree_, position_);
      // run computation for the i-th coefficient
      if (!coefs_[i].run())
      {
        msg_error_ =coefs_[i].error();
        return false;
      }
      // get coefficients
      coefficients_.merge(coefs_[i].coefficients());
    }
  }
  catch ( Exception const& e)
  {
    msg_error_ = e.error();
    return false;
  }
#ifdef STK_REGRESS_VERBOSE
  stk_cout << _T("AdditiveBSplineCoefficients::run() done\n");
#endif
  return true;
}

/* run the computations. */
template<class Array>
ArrayXX AdditiveBSplineCoefficients<Array>::extrapolate(Array const& x) const
{
#ifdef STK_REGRESS_VERBOSE
  stk_cout << _T("in AdditiveBSplineCoefficients::extrapolate()\n");
#endif
  if (x.cols() != coefs_.range())
  { STKRUNTIME_ERROR_NO_ARG(AdditiveBSplineBSplineCoefficients::extrapolate(x),x.cols() != coefs_.range());}
    // resize the matrix of coefficient
    ArrayXX coefficients(x.rows(), Range());
    for (int j= coefs_.begin(); j<coefs_.end(); j++)
    { coefficients.merge(coefs_[j].extrapolate(x.col(j)));}
#ifdef STK_REGRESS_VERBOSE
  stk_cout << _T("AdditiveBSplineCoefficients::run() done\n");
#endif
  return coefficients;
}


/* Compute the coefficients of the B-spline curve for the given values.
 *  @param p_data the input data values
 *  @param nbControlPoints number of control points
 *  @param degree degree of the B-spline curves
 *  @param position method to use for positioning the knots
 **/
template<class Array>
void AdditiveBSplineCoefficients<Array>::setData( Array const* p_data
                                               , int nbControlPoints
                                               , int degree
                                               , Regress::KnotsPosition position
                                               )
{
  p_data_ =p_data;
  nbKnots_ = nbControlPoints + degree +1;
  nbControlPoints_ = nbControlPoints;
  degree_ = degree;
  position_ = position;
}


} // namespace STK

#endif /* STK_BSPLINECOEFFICIENTS_H */
