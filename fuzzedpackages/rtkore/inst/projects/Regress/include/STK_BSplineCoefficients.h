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
 * Purpose:  Compute the coefficient of a B-spline curves.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_BSplineCoefficients.h
 *  @brief In this file we define the BSplineCoefficients class.
 **/

#ifndef STK_BSPLINECOEFFICIENTS_H
#define STK_BSPLINECOEFFICIENTS_H

#include "STK_IBasis.h"

#include <Arrays/include/STK_Array2DVector.h>
#include <DManager/include/STK_HeapSort.h>

#ifdef STK_REGRESS_VERBOSE
#include <Arrays/include/STK_Display.h>
#endif

namespace STK
{
/** @brief Compute the regression B-splines coefficients.
 * The BSplineCoefficients class computes the coefficients of a BSpline curve
 * using the de Boor's algorithm. The knots can be uniform (the default),
 * periodic or density placed.
 *
 * Given the number of control points and the degree of the B-splines, the number
 * of knots is given by the formula @e nbKnots = nbControlPoints + degree + 1.
 *
 * @note If the input data set is a vector of size @c n the output matrix of the
 * coefficients @c Coefficients() is a matrix of size @c (n, nbControlPoints).
 */
template<class Data>
class BSplineCoefficients: public IBasis<Data, ArrayXX>
{
  public:
    typedef IBasis<Data, ArrayXX> Base;
    using Base::p_data_;
    using Base::coefficients_;
    using Base::dim_;
    using Base::msg_error_;
    using Base::hasRun_;
    using Base::minValue_;
    using Base::maxValue_;

    /** @brief Default constructor : initialize the data members with default
     *  values.
     *  The number of knots is given by the formula
     *  nbKnots = nbControlPoints + degree +1.
     *  @param p_data the input data values
     *  @param nbControlPoints number of control points
     *  @param degree degree of the B-spline curves
     *  @param position method to use for positioning the knots
     *  @param useDataValues if @c true use data in order to find the minValue_ and maxValue_, use
     *  values set otherwise.
     **/
    BSplineCoefficients( Data const* p_data =0
                       , int nbControlPoints =1
                       , int degree = 3
                       , Regress::KnotsPosition position = Regress::uniformKnotsPositions_
                       , bool useDataValues = true
                       );
    /** Constructor : initialize the data members. The number of knots is given
     *  by the formula nbKnots = nbControlPoints + degree +1.
     *  @param data the input data values
     *  @param nbControlPoints number of control points
     *  @param degree degree of the B-spline curves
     *  @param position method to use for positioning the knots
     *  @param useDataValues if @c true use data in order to find the minValue_ and maxValue_, use
     *  values set otherwise.
     **/
    BSplineCoefficients( Data const& data
                       , int nbControlPoints
                       , int degree = 3
                       , Regress::KnotsPosition position = Regress::uniformKnotsPositions_
                       , bool useDataValues = true
                       );
    /** copy constructor.
     *  @param coefs the coefficients to copy
     **/
    BSplineCoefficients( BSplineCoefficients const& coefs);
    /** Destructor. */
    virtual ~BSplineCoefficients() {}
    /** clone pattern implementation */
    BSplineCoefficients* clone() const { return new BSplineCoefficients(*this);}
    /** run the computations. */
    virtual bool run();
    // getters
    /** @return the degree of the B-Splines */
    inline int degree() const { return degree_;}
    /** @return the number of knots of the B-spline curve */
    inline int nbKnots() const { return nbKnots_;}
    /** @return the number of control points of the curve (as dim())*/
    inline int nbControlPoints() const { return dim_;}
    /** @return the vector of knots of the B-spline curve */
    inline VectorX const& knots() const { return knots_;}
    // setters
    /** Set the number of control point (the number of BSpline)
     *  @param nbControlPoints number of control points
     **/
    void setNbControlPoints( int nbControlPoints);
    /** Alias for setting the number of control point (the number of BSpline)
     *  @param dim number of control points
     **/
    inline void setDim( int dim) { setNbControlPoints(dim);}
    /** Set the degree of the BSpline basis
     *  @param degree degree of the B-spline curves (default is 3)
     **/
    void setDegree( int degree = 3);
    /** Set the kind of position to use for the knots
     *  @param position method to use for positioning the knots (default is uniform)
     **/
    void setPosition(Regress::KnotsPosition position);
    /**Set the whole set of parameters for computing the coefficients
     *  @param data the input data values
     *  @param nbControlPoints number of control points
     *  @param degree degree of the B-spline curves (default is 3)
     *  @param position method to use for positioning the knots (default is uniform)
     **/
    void setParameters( Data const& data
                      , int nbControlPoints
                      , int degree = 3
                      , Regress::KnotsPosition position = Regress::uniformKnotsPositions_
                      );

    /** [DEPRECATED] Set the whole set of parameters for computing the coefficients
     *  @param data the input data values
     *  @param nbControlPoints number of control points
     *  @param degree degree of the B-spline curves (default is 3)
     *  @param position method to use for positioning the knots (default is uniform)
     **/
    void setData( Data const& data
                , int nbControlPoints
                , int degree = 3
                , Regress::KnotsPosition position = Regress::uniformKnotsPositions_
                );

    /** @return the extrapolated matrix of coefficients for a given set of x-values.
     *  @param x the values to extrapolate
     **/
    template<class OtherVector>
    ArrayXX extrapolate(OtherVector const& x) const;

  protected:
    /** degree of the B-splines curves. */
    int degree_;
    /** Method used in order to position the knots. */
    Regress::KnotsPosition position_;
    /** number of knots of the B-spline curves.*/
    int nbKnots_;
    /** Index of the last control point of the B-spline curves.
     *  This is dim_ - 1.
     **/
    int lastControlPoint_;
    /** Data of the knots */
    VectorX knots_;

    /** compute the position of the knots of the B-spline curves.*/
    bool computeKnots();
    /** Compute the coefficients of the B-spline curves.*/
    void computeCoefficients();

  private:
    /** compute the position of the uniform knots.*/
    void computeUniformKnots();
    /** compute the position of the periodic knots.*/
    void computePeriodicKnots();
    /** compute the position of the density knots
     *  @param isSorted @c true if the data is sorted, @c false otherwise
     **/
    void computeDensityKnots(bool isSorted);
    /** Compute a row of the coefficients matrix for a given value.
     * @param irow index of the row
     * @param value the value to which we want to compute the row
     **/
    void computeCoefficientsRow(int irow, Real const& value);
};

/* constructor */
template<class Data>
BSplineCoefficients<Data>::BSplineCoefficients( Data const* p_data
                                                , int nbControlPoints
                                                , int degree
                                                , Regress::KnotsPosition position
                                                , bool useDataValues
                                                )
                                                : Base(p_data, nbControlPoints, useDataValues)
                                                , degree_(degree)
                                                , position_(position)
                                                , nbKnots_(nbControlPoints + degree + 1)
                                                , lastControlPoint_(dim_-1)
                                                , knots_( Range(0, nbKnots_) )
{ }

/* constructor */
template<class Data>
BSplineCoefficients<Data>::BSplineCoefficients( Data const& data
                                              , int nbControlPoints
                                              , int degree
                                              , Regress::KnotsPosition position
                                              , bool useDataValues
                                              )
                                              : Base(data, nbControlPoints, useDataValues)
                                              , degree_(degree)
                                              , position_(position)
                                              , nbKnots_(nbControlPoints + degree +1)
                                              , lastControlPoint_(dim_-1)
                                              , knots_( Range(0, nbKnots_) )
{}
/* copy constructor.
 *  @param coefs the coefficients to copy
 **/
template<class Data>
BSplineCoefficients<Data>::BSplineCoefficients( BSplineCoefficients const& coefs)
                                              : Base(coefs)
                                              , degree_(coefs.degree_)
                                              , position_(coefs.position_)
                                              , nbKnots_(coefs.nbKnots_)
                                              , lastControlPoint_(coefs.lastControlPoint_)
                                              , knots_(coefs.knots_)
{}



/*  run the computations for the given value.
 *  @param p_data the input data values
 **/
template<class Data>
void BSplineCoefficients<Data>::setData( Data const& data
                                       , int nbControlPoints
                                       , int degree
                                       , Regress::KnotsPosition position
                                       )
{ // set data
  p_data_ = &data;
  dim_ = nbControlPoints;
  degree_ = degree;
  position_ = position;
  nbKnots_ = dim_ + degree_ +1;
  lastControlPoint_ = dim_-1;
  Base::update();
}

/*  run the computations for the given value.
 *  @param p_data the input data values
 **/
template<class Data>
void BSplineCoefficients<Data>::setParameters( Data const& data
                                             , int nbControlPoints
                                             , int degree
                                             , Regress::KnotsPosition position
                                             )
{ // set data
  p_data_ = &data;
  dim_ = nbControlPoints;
  degree_ = degree;
  position_ = position;
  nbKnots_ = dim_ + degree_ +1;
  lastControlPoint_ = dim_-1;
  knots_.resize( Range(0, nbKnots_) );
  Base::update();
}

/*  run the computations for the given value.
 *  @param p_data the input data values
 **/
template<class Data>
void BSplineCoefficients<Data>::setNbControlPoints( int nbControlPoints)
{
  dim_ = nbControlPoints;
  nbKnots_ = dim_ + degree_ +1;
  lastControlPoint_ = dim_-1;
  knots_.resize( Range(0, nbKnots_) );
  Base::update();
}

/* Set the degree of the BSpline basis
 *  @param degree degree of the B-spline curves (default is 3)
 **/
template<class Data>
void BSplineCoefficients<Data>::setDegree( int degree)
{
  degree_ = degree;
  nbKnots_ = dim_ + degree_ +1;
  lastControlPoint_ = dim_-1;
  knots_.resize( Range(0, nbKnots_)  );
  Base::update();
}

/* Set the kind of position to use for the knots
 *  @param position method to use for positioning the knots (default is uniform)
 **/
template<class Data>
void BSplineCoefficients<Data>::setPosition(Regress::KnotsPosition position)
{
  position_ = position;
  Base::update();
}


/* run the computations. */
template<class Data>
bool BSplineCoefficients<Data>::run()
{
  // check if data exists
  if (!p_data_)
  {
   msg_error_ = STKERROR_NO_ARG(Error in BSplineCoefficients::run,p_data_ is null);
   return false;
  }
  if (!this->initializeStep()) return false;
  knots_ = minValue_;
  // compute the knots and coefficients
  if (computeKnots()) { computeCoefficients();}
  else                { return false;}
  this->hasRun_ = true;
  return true;
}

/* Extrapolate the matrix of coefficients for a given set of x-values.
 *  @param x the values to extrapolate
 *  @param coefs the matrix of coefficients for each values.
 **/
template<class Data>
template<class OtherVector>
ArrayXX BSplineCoefficients<Data>::extrapolate(OtherVector const& x) const
{
  // check if knots exists
  if (!this->hasRun_)
  { STKRUNTIME_ERROR_NO_ARG(BSplineCoefficients::extrapolate,No run);}
  // resize coeficients
  ArrayXX coefs(x.range(), Range(0, lastControlPoint_, 0), 0.0);
  // check if the original data set was not reduced to a single point
  if (minValue_ == maxValue_) return coefs;
  // compute the coefficients
  for (int irow=x.begin(); irow< x.end(); irow++)
  {
    const Real value = x[irow];
    // value outside the range of the knots case
    if (value <= minValue_)
    {
      coefs(irow, 0) = 1.0;
      continue;
    }
    if (value >= maxValue_)
    {
      coefs(irow, lastControlPoint_) = 1.0;
      continue;
    }
    // find interval
    int k, k1;
    for (k=0, k1=1; k<lastControlPoint_; k++, k1++)
    {
      if (value < knots_[k1]) break;
    }
    // begin recursion
    coefs(irow, k) = 1.0;
    for (int d=1; d<=degree_; d++)
    {
      // right (south-west corner) term only
      coefs(irow, k-d) = ( (knots_[k1] - value)/(knots_[k1] - knots_[k1-d]) ) * coefs(irow, k1-d);
      // compute internal terms
      for (int i = k1-d; i<k; i++)
      {
        const Real knots_i = knots_[i], knots_id1 = knots_[i+d+1];
        coefs(irow, i) = ( (value - knots_i)/(knots_[i+d] - knots_i) ) * coefs(irow, i)
                       + ( (knots_id1 - value)/(knots_id1 - knots_[i+1]) ) * coefs(irow, i+1);
      }
      // left (north-west corner) term only
      coefs(irow, k) *= (value - knots_[k])/(knots_[k+d] - knots_[k]);
    }
  }
  return coefs;
}

/* compute the knots of the B-spline curves.*/
template<class Data>
bool BSplineCoefficients<Data>::computeKnots()
{
  // resize and initialize knots
  knots_.resize( Range(0, nbKnots_) ) = minValue_;
  // set knots values
  switch (position_)
  {
    // uniform position
    case Regress::uniformKnotsPositions_:
      computeUniformKnots();
      break;
    // periodic position
    case Regress::periodicKnotsPositions_:
      computePeriodicKnots();
      break;
    // density position
    case Regress::densityKnotsPositions_:
      computeDensityKnots(true);
      break;
    default:

      msg_error_ = STKERROR_NO_ARG(BSplineCoefficients::computeKnots,invalid position);
      return false;
      break;
  }
  // shift knots
  if (position_ != Regress::densityKnotsPositions_)
  {
    Real range = (maxValue_ - minValue_);
    for (int k = 0; k < nbKnots_; k++) knots_[k] = minValue_ + range * knots_[k];
  }
  return true;
}

/* Compute the coefficients of the B-spline curves.*/
template<class Data>
void BSplineCoefficients<Data>::computeCoefficients()
{
#ifdef STK_REGRESS_VERBOSE
  stk_cout << _T("BSplineCoefficients::computeCoefficients()\n");
#endif

  // compute the coefficients
  for (int i=p_data_->begin(); i< p_data_->end(); i++)
  { computeCoefficientsRow(i, (*p_data_)[i]);}

#ifdef STK_REGRESS_VERBOSE
  stk_cout << _T("BSplineCoefficients::computeCoefficients() done\n");
#endif
}

/* compute the position of the uniform knots.*/
template<class Data>
void BSplineCoefficients<Data>::computeUniformKnots()
{
  // compute step
  Real step = 1.0/(dim_ - degree_);
  // set internal knots
  const int first = degree_ + 1;
  for (int k = first, j = 1; k <= lastControlPoint_; j++, k++)  knots_[k] = j * step;
  // set external knots
  for ( int k=0, j = lastControlPoint_+1; k < first; j++, k++)
  {
    knots_[k] = 0;
    knots_[j] = 1;
  }
}
/* compute the position of the periodic knots.*/
template<class Data>
void BSplineCoefficients<Data>::computePeriodicKnots()
{
  // compute step
  Real step = 1.0/(dim_ - degree_);
  // set knots
  for (int k = 0, j = -degree_; k < nbKnots_; j++, k++)
    knots_[k] = j * step;
;
}
/* compute the position of the density knots. */
template<class Data>
void BSplineCoefficients<Data>::computeDensityKnots(bool isSorted)
{
  // sorted data
  Data xtri(*p_data_, true);
  // sort the data
  if (!isSorted) heapSort< Data >(xtri);

  // compute step
  Real step = xtri.size()/(Real)(lastControlPoint_-degree_+1);

  // set internal knots
  int first = xtri.begin();
  for (int k = degree_ + 1, kcell =1; k <= lastControlPoint_; k++, kcell++)
  {
    int cell = first + int(kcell* step);
    knots_[k] = (xtri[cell] + xtri[cell+1])/2.;
  }
  // set external knots
  for ( int k=0, j = lastControlPoint_+1; k <= degree_; j++, k++)
  {
    knots_[k] = xtri.front();
    knots_[j] = xtri.back();
  }
}

/* Compute a row of the coefficients
 * @param irow index of the row
 **/
template<class Data>
void BSplineCoefficients<Data>::computeCoefficientsRow(int irow, Real const& value)
{
  // value outside the range of the knots case
  if (value <= minValue_)
  {
    coefficients_(irow, 0) = 1.0;
    return;
  }
  if (value >= maxValue_)
  {
    coefficients_(irow, lastControlPoint_) = 1.0;
    return;
  }
  // find interval
  int k, k1;
  for (k=0, k1=1; k<lastControlPoint_; k++, k1++)
  {
    if (value < knots_[k1]) break;
  }
  // begin recursion
  coefficients_(irow, k) = 1.0;
  for (int d=1; d<=degree_; d++)
  {
    // right (south-west corner) term only
    coefficients_(irow, k-d) = ( (knots_[k1] - value)/(knots_[k1] - knots_[k1-d]) )
                               * coefficients_(irow, k1-d);
    // compute internal terms
    for (int i = k1-d; i<k; i++)
    {
      const Real knots_i = knots_[i], knots_id1 = knots_[i+d+1];
      coefficients_(irow, i) = ( (value - knots_i)/(knots_[i+d] - knots_i) )
                               * coefficients_(irow, i)
                             + ( (knots_id1 - value)/(knots_id1 - knots_[i+1]) )
                               * coefficients_(irow, i+1);
    }
    // left (north-west corner) term only
    coefficients_(irow, k) *= (value - knots_[k])/(knots_[k+d] - knots_[k]);
  }
}

} // namespace STK

#endif /* STK_BSPLINECOEFFICIENTS_H */
