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

/*    Project: stkpp::Regress
 * created on: Nov 28, 2017
 *     Author: iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IBasis.h
 *  @brief In this file we define the Interface class IBasis for basis functions
 **/

#ifndef STK_IBASIS_H
#define STK_IBASIS_H

#include "STK_Regress_Util.h"

#include <Sdk/include/STK_IRunner.h>
#include <Arrays/include/STK_Array2D.h>

namespace STK
{

/** @ingroup Regress
 *  Interface base class for all basis function
 **/
template<class Data, class Coefs = ArrayXX>
class IBasis: public IRunnerBase
{
  public:
    typedef typename Data::Type Type;
    /** default constructor */
    IBasis();
    /** constructor
     *  @param p_data pointer on the data set (can be null)
     *  @param dim dimension of the basis
     *  @param useDataValues if @c true use data in order to find the minValue_ and maxValue_, use
     *  values set otherwise.
     **/
    IBasis(Data const* p_data, int dim, bool useDataValues = true);
    /** constructor
     *  @param data reference on the data set
     *  @param dim dimension of the basis
     *  @param useDataValues if @c true use data in order to find the minValue_ and maxValue_, use
     *  values set otherwise.
     **/
    IBasis(Data const& data, int dim, bool useDataValues = true);
    /** copy constructor.
     *  @param basis the basis to copy
     **/
    IBasis( IBasis const& basis);
    /** destructor */
    virtual ~IBasis() {}

    // getters
    /** @return the dimension of the basis */
    inline int dim() const { return dim_;}
    /** @return the matrix with the coefficients of the basis curve. */
    inline Coefs const& coefficients() const { return coefficients_;}
    /** @return the minimal value of the data */
    inline Type minValue() const { return minValue_;}
    /** @return the maximal value of the data */
    inline Type maxValue() const { return maxValue_;}

    // setters
    /** Set the data set
     *  @param data the input data values
     **/
    inline void setData( Data const& data)
    {
      p_data_ = &data;
      update();
    }
    /** @param dim number of dimension */
    inline void setDim( int dim) { dim_ = dim; update();}
    /** @param minValue minimal value of the data */
    inline void setMinValue( Type const& minValue) { minValue_ = minValue;}
    /** @param maxValue maximal value of the data */
    inline void setMaxValue( Type const& maxValue) { maxValue_ = maxValue;}

    /** Initialize the parameters minValue_ and maxValue_ using data set */
    bool initializeStep();

  protected:
    /** update IBasis
     *  if a parameter or a new data set is set, update the state of this runner.
     **/
    virtual void update();
    /** the input data set */
    Data const* p_data_;
    /** number of dimension to build*/
    int dim_;
    /** */
    bool useDataValues_;
    /** Minimal value of the data */
    Type minValue_;
    /** Maximal value of the data */
    Type maxValue_;

    /** Array2D of the coefficients */
    Coefs coefficients_;
};

template<class Data, class Coefs>
IBasis<Data, Coefs>::IBasis(): IRunnerBase()
                            , p_data_(0)
                            , dim_()
                            , useDataValues_(true)
                            , minValue_( Arithmetic<Type>::max())
                            , maxValue_(-Arithmetic<Type>::max())
                            , coefficients_()
{}

template<class Data, class Coefs>
IBasis<Data, Coefs>::IBasis( Data const* p_data, int dim, bool useDataValues)
                           : IRunnerBase()
                           , p_data_(p_data)
                           , dim_(dim)
                           , useDataValues_(useDataValues)
                           , minValue_( Arithmetic<Type>::max())
                           , maxValue_(-Arithmetic<Type>::max())
                           , coefficients_()
{}
template<class Data, class Coefs>
IBasis<Data, Coefs>::IBasis( Data const& data, int dim, bool useDataValues)
                           : IRunnerBase()
                           , p_data_(&data)
                           , dim_(dim)
                           , useDataValues_(useDataValues)
                           , minValue_( Arithmetic<Type>::max())
                           , maxValue_(-Arithmetic<Type>::max())
                           , coefficients_()
{}

template<class Data, class Coefs>
IBasis<Data, Coefs>::IBasis( IBasis const& basis)
                           : IRunnerBase(basis)
                           , p_data_(basis.p_data_)
                           , dim_(basis.dim_)
                           , useDataValues_(basis.useDataValues_)
                           , minValue_(basis.minValue_)
                           , maxValue_(basis.maxValue_)
                           , coefficients_(basis.coefficients_)
{}

/* update IBasis
 * if a parameter or a new data set is set, update the values of this runner.
 **/
template<class Data, class Coefs>
void IBasis<Data, Coefs>::update()
{
  if (useDataValues_)
  {
    minValue_ =  Arithmetic<Type>::max();
    maxValue_ = -Arithmetic<Type>::max();
  }
  hasRun_ = false;
}
/* Initialize the parameters */
template<class Data, class Coefs>
bool IBasis<Data, Coefs>::initializeStep()
{
  // resize and initialize coeficients
  coefficients_.resize(p_data_->range(), Range(0, dim_)) =0;
  // compute min and max value
  if (useDataValues_)
  {
    minValue_ =  Arithmetic<Type>::max();
    maxValue_ = -Arithmetic<Type>::max();
    for (int i=p_data_->begin(); i< p_data_->end(); i++)
    {
      minValue_ = std::min(minValue_, (*p_data_)[i]);
      maxValue_ = std::max(maxValue_, (*p_data_)[i]);
    }
  }
  // if all value are equals
  if (minValue_ == maxValue_)
  {
    msg_error_ = STKERROR_NO_ARG(IBasis::initializeStep,All values are equal);
    return false;
  }
  // if values are incorrect
  if (minValue_ > maxValue_)
  {
    msg_error_ = STKERROR_NO_ARG(IBasis::initializeStep,minValue_ greater than maxValue_);
    return false;
  }
  return true;
}


} // namespace STK

#endif /* STK_IBASIS_H_ */
