/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2017  Serge Iovleff, Université Lille 1, Inria

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
 * created on: Oct 26, 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ChebyshevCoefficients.h
 *  @brief In this file we define the ChebyshevCoefficients class
 **/

#ifndef STK_CHEBYSHEVCOEFFICIENTS_H
#define STK_CHEBYSHEVCOEFFICIENTS_H

#include "STK_IBasis.h"

namespace STK
{

/** @brief ChebyshevCoefficients class allows to compute the coefficients of a sampled function
 *  using Chebyshev polynomials.
 *
 *  Chebyshev polynomials are important in approximation theory because the roots of the Chebyshev
 *  polynomials of the first kind, which are also called Chebyshev nodes, are used as nodes in
 *  polynomial interpolation.
 */
template<class Data, class Coefs = ArrayXX>
class ChebyshevCoefficients: public IBasis<Data, Coefs>
{
  public:
    typedef IBasis<Data, Coefs> Base;
    typedef typename Data::Type Type;

    using Base::p_data_;
    using Base::coefficients_;
    using Base::dim_;
    using Base::minValue_;
    using Base::maxValue_;
    using Base::msg_error_;
    using Base::hasRun_;

    /** @brief default constructor
     *  @param p_data pointer on the data set
     *  @param dim number of basis to use
     *  @param useDataValues if @c true use data in order to find the minValue_ and maxValue_, use
     *  values set otherwise
     **/
    ChebyshevCoefficients( Data const* p_data =0, int dim =1, bool useDataValues = true)
                         : Base(p_data, dim, useDataValues){}
    /** @brief constructor
     *  @param data reference on the data set
     *  @param dim number of basis to use
     *  @param useDataValues  if @c true use data in order to find the minValue_ and maxValue_, use
     *  values set otherwise
     **/
    ChebyshevCoefficients( Data const& data, int dim, bool useDataValues = true)
                         : Base(data, dim, useDataValues){}
    /** copy constructor.
     *  @param coefs the coefficients to copy
     **/
    ChebyshevCoefficients( ChebyshevCoefficients const& coefs) : Base(coefs){}
    /** Destructor. */
    virtual ~ChebyshevCoefficients() {}
    /** clone pattern implementation */
    inline ChebyshevCoefficients* clone() const { return new ChebyshevCoefficients(*this);}

    /** run the computations. */
    virtual bool run();
};

template<class Data, class Coefs>
bool ChebyshevCoefficients<Data,Coefs>::run()
{
  // check if data exists
  if (!p_data_)
  {
   msg_error_ = STKERROR_NO_ARG(Error in SinesCoefficients::run,p_data_ is null);
   return false;
  }
  if (!this->initializeStep()) return false;
  Data x = ( Type(2)* (*p_data_)-(minValue_+maxValue_))/(maxValue_-minValue_);
  // resize and initialize coefficients
  coefficients_.resize(p_data_->range(), Range(0, dim_)) =0;
  if (dim_>1)
  {
    coefficients_.col(0) = Type(1);
    if (dim_>2)
    {
      coefficients_.col(1) = x;
    }
    // compute T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x)
    for(int j=2; j<coefficients_.endCols(); ++j)
    {
      coefficients_.col(j) = 2. * x *coefficients_.col(j-1) - coefficients_.col(j-2); ;
    }
  }
  // rescale
  for(int d=1; d<coefficients_.endCols(); ++d)
  {
    coefficients_.col(d) = (coefficients_.col(d) * (maxValue_-minValue_) + (minValue_+maxValue_))/Type(2);
  }
  this->hasRun_ = true;
  return true;

}

} // namespace STK

#endif /* STK_CHEBYSHEVCOEFFICIENTS_H */
