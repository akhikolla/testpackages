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

/** @file STK_SinesCoefficients.h
 *  @brief In this file we define the SinesCoefficients class
 **/

#ifndef STK_SINUSCOEFFICIENTS_H
#define STK_SINUSCOEFFICIENTS_H

#include "STK_IBasis.h"

namespace STK
{

/** @brief SinesCoefficients class allows to compute the coefficients of a sampled function
 *  using sines basis functions
 */
template<class Data, class Coefs = ArrayXX>
class SinesCoefficients: public IBasis<Data, Coefs>
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
     *  values set otherwise.
     **/
    SinesCoefficients( Data const* p_data =0, int dim =1, bool useDataValues = true)
                     : Base(p_data, dim, useDataValues){}
    /** @brief constructor
     *  @param data reference on the data set
     *  @param dim number of basis to use
     *  @param useDataValues if @c true use data in order to find the minValue_ and maxValue_, use
     *  values set otherwise.
     **/
    SinesCoefficients( Data const& data, int dim, bool useDataValues = true)
                     : Base(data, dim, useDataValues){}
    /** copy constructor.
     *  @param coefs the coefficients to copy
     **/
    SinesCoefficients( SinesCoefficients const& coefs) : Base(coefs){}
    /** Destructor. */
    virtual ~SinesCoefficients() {}
    /** clone pattern implementation */
    inline SinesCoefficients* clone() const { return new SinesCoefficients(*this);}

    /** run the computations. */
    virtual bool run();
};

template<class Data, class Coefs>
bool SinesCoefficients<Data,Coefs>::run()
{
  // check if data exists
  if (!p_data_)
  {
   msg_error_ = STKERROR_NO_ARG(Error in SinesCoefficients::run,p_data_ is null);
   return false;
  }
  if (!this->initializeStep()) return false;
  if (dim_>1)
  {
    coefficients_.col(0) = Type(1);
    Type period = Type(2)*Const::_PI_/(maxValue_ - minValue_);
    for(int j=1; j<coefficients_.endCols(); ++j)
    {
      coefficients_.col(j) = (*p_data_ * Type(j)*period).sin();
    }
  }
  this->hasRun_ = true;
  return true;

}

} // namespace STK

#endif /* STK_SINUSCOEFFICIENTS_H */
