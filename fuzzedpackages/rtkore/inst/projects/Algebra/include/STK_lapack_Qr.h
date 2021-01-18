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
 * Project:  stkpp::Algebra
 * created on: 20 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_lapack_Qr.h
 *  @brief In this file we define the enclosing class of the dgeqrf lapack routine.
 **/


#ifndef STK_LAPACK_QR_H
#define STK_LAPACK_QR_H

#include <Arrays/include/STK_CArray.h>
#include <Arrays/include/STK_CArrayVector.h>
#include <Arrays/include/STK_Array2D.h>

#include "STK_IQr.h"
#include "STK_lapack_Util.h"

namespace STK
{

namespace lapack
{
class Qr;
}

namespace hidden
{
/** @ingroup hidden
 *  Specialization for the Qr class.
 **
 **/
template<>
struct AlgebraTraits< lapack::Qr >
{
  typedef ArrayXX Array;
  typedef typename Array::Col ColVector;
  typedef typename Array::Row RowVector;
};

} // namespace hidden


namespace lapack
{
/** @ingroup Algebra
 *  {
 *    @class Qr
 *    @brief Qr computes the QR decomposition of a real matrix using the
 *    Lapack routine *geqrf.
 *
 *  - Input:  A matrix (nrow,ncol)
 *  - Output:
 *    -# Q  matrix (nrow,ncol) with the Housholder vectors in the min(nrow, ncol) first columns.
 *    -# R  matrix (nrow,ncol) upper triangular.
 *
 *    @sa STK::IQr
 *
 */
class Qr: public IQr<Qr >
{
  public:
    typedef IQr<Qr > Base;
    using Base::Q_;
    using Base::R_;
    /** Default constructor.
     *  @param data the matrix to decompose
     *  @param ref true if we overwrite A
     **/
    inline Qr( ArrayXX const&  data, bool ref = false): Base(data, ref) {}
    /** @brief Constructor
     *  @param data reference on a matrix expression
     */
    template<class Derived>
    Qr( ArrayBase<Derived> const& data): Base(data){}
    /** Copy constructor.
     *  @param decomp the decomposition  to copy
     **/
    inline Qr( Qr const& decomp): Base(decomp) {}
    /** virtual destructor */
    inline virtual ~Qr() {}
    /** @brief clone pattern */
    inline virtual Qr* clone() const { return new Qr(*this);}
    /** Operator = : overwrite the Qr with decomp. */
    inline Qr& operator=(Qr const& decomp)
    {
      Base::operator =(decomp);
      return *this;
    }
    /** @brief Run qr decomposition
     *  Launch geqrf LAPACK routine to perform the qr decomposition.
     *  @return @c true if no error occur, @c false otherwise
     */
    bool runImpl();

  private:
    /** private method for computing the Qr decomposition using a CArrayXX array */
    bool computeQr(CArrayXX& a, CVectorX& tau);
};
/** @} */

/* private method for computing the Qr decomposition using a CArrayXX array */
inline bool Qr::computeQr(CArrayXX& a, CVectorX& tau)
{
  // start qr iterations
  int lWork =-1, m = a.sizeRows(), n= a.sizeCols();
  int info;
  Real iwork;
  Real *p_work;

  // get size for work
  info = geqrf(m, n, a.p_data(), m, tau.p_data(), &iwork, lWork);
  // check any error
  if (info!=0)
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(lapack::Qr::computeQr,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(lapack::Qr::computeQr get,-info,error parameter);
    return false;
  }
  // create
  lWork = (int)iwork;
  p_work = new Real[lWork];
  info = geqrf(m, n, a.p_data(), m, tau.p_data(), p_work, lWork);
  // release working set
  delete[] p_work;
  // check any error
  if (info!=0)
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(lapack::Qr::computeQr,internal error);
     return false;
    }
    this->msg_error_= STKERROR_1ARG(lapack::Qr::computeQr get,-info,error parameter);
    return false;
  }
  return true;
}
/* Computing the QR decomposition of the matrix Q_. */
inline bool Qr::runImpl()
{
  int begin = Q_.beginRows(); // same first index checked at construction of QR
  int end   = std::min(Q_.endRows(), Q_.endCols());
  // compute results
  CArrayXX a = Q_;
  CVectorX tau(std::min(Q_.sizeRows(), Q_.sizeCols()));
  a.shift(0, 0);
  tau.shift(0);
  if (!computeQr(a, tau)) { return false;};
  // get results
  a.shift(begin, begin);
  tau.shift(begin);
  R_.resize(Q_.rows(), Q_.cols());
  for (int i=Q_.beginRows(); i< Q_.endRows(); ++i)
  {
    int jmax = std::min(i, Q_.endCols());
    for (int j=Q_.beginCols(); j< jmax; ++j) { Q_(i,j) = a(i,j);}
    for (int j=i; j< Q_.endCols(); ++j)      { R_(i,j) = a(i,j);}
  }
  for (int j=Q_.beginCols(); j<end; ++j) { Q_(j,j) = -tau[j];}
  return true;
}


} // namespace lapack

} // namespace STK

#endif /* STK_LAPACK_QR_H */
