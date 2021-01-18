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
 * Purpose:  Define the Qr Class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Qr.h
 *  @brief In thThis file wedefine and implement the Qr class.
 **/

#ifndef STK_QR_H
#define STK_QR_H

#include "STK_IQr.h"
#include "STK_Givens.h"
#include <Arrays/include/STK_Array2D.h>
#include <Arrays/include/STK_Array2DPoint.h>

#ifdef STK_ALGEBRA_DEBUG
#include <Arrays/include/STK_Display.h>

template< class Container2D >
void print(Container2D const& A, STK::String const& name)
{
  stk_cout << "print: " << name << _T("\n";);
  stk_cout << name << _T(".isRef() =")        << A.isRef()  << _T("\n");
  stk_cout << name << _T(".cols() =")      << A.cols()  << _T("\n");
  stk_cout << name << _T(".rows() =")      << A.rows()  << _T("\n\n");
}

void print(ArrayXX const& A, STK::String const& name)
{
  stk_cout << "print: " << name << _T("\n";);
  stk_cout << name << _T(".isRef() =")        << A.isRef()  << _T("\n");
  stk_cout << name << _T(".cols() =")      << A.cols()  << _T("\n");
  stk_cout << name << _T(".rows() =")      << A.rows()  << _T("\n\n");
  stk_cout << name << _T(".rangeCols().isRef() =")  << A.rangeCols().isRef() << _T("\n");
  stk_cout << name << _T(".rangeCols() =\n")  << A.rangeCols() << _T("\n");
}
#endif

namespace STK
{
// forward declaration
class Qr;

namespace hidden
{
/** @ingroup hidden
 *  Specialization for the Qr class.
 **/
template<>
struct AlgebraTraits< Qr >
{
  typedef ArrayXX Array;
  typedef typename Array::Col ColVector;
  typedef typename Array::Row RowVector;
};

} // namespace hidden

/** @ingroup Algebra
 *  @brief The class Qr perform the QR decomposition of an ArrayXX.
 *
 *  - Input:  A matrix (nrow,ncol)
 *  - Output:
 *    -# Q  matrix (nrow,ncol) with the Housholder vectors in the min(nrow, ncol) first columns.
 *    -# R  matrix (nrow,ncol) upper triangular.
 *
 *  @sa STK::IQr
 */
class Qr: public IQr<Qr >
{
  public :
    typedef typename hidden::AlgebraTraits< Qr >::ColVector ColVector;
    typedef typename hidden::AlgebraTraits< Qr >::RowVector RowVector;
    typedef IQr<Qr> Base;
    using Base::Q_;
    using Base::R_;

    /** Default constructor.
     *  @param A the matrix to decompose
     *  @param ref true if we overwrite A
     **/
    inline Qr( ArrayXX const& A, bool ref = false): Base(A, ref) {}
    /** @brief Constructor
     *  @param data reference on a matrix expression
     */
    template<class Derived>
    inline Qr( ExprBase<Derived> const& data): Base(data) {}
    /** Copy constructor.
     *  @param decomp the decomposition  to copy
     **/
    inline Qr( Qr const& decomp): Base(decomp) {}
    /** virtual destructor */
    inline virtual ~Qr() {}
    /** clone pattern */
    inline virtual Qr* clone() const { return new Qr(*this);}
    /** Operator = : overwrite the Qr with decomp. */
    inline Qr& operator=(Qr const& decomp)
    { Base::operator=(decomp);
      return *this;
    }
    /** Compute the QR decomposition. **/
    inline bool runImpl() { qr(); return true;}

  private:
  /** Compute the qr decomposition of the matrix Q_ */
    void qr();
};

/* Computation of the QR decomposition */
inline void Qr::qr()
{
  R_.resize(Q_.rows(), Q_.cols());
  // start left householder reflections
  Range r(Q_.rows()), c(Q_.cols());
  for(int j = R_.beginRows(); (j < R_.endRows()) && (j < R_.endCols()) ; ++j)
  {
    ColVector u(Q_, r, j);// get a reference on the j-th column in the range r
    R_(j, j) = house(u);  // compute the housolder vector
    applyLeftHouseholderVector(Q_.sub(r, c.incFirst(1)), u); // apply-it to the remaining part of Q_
    R_.row(j, c).copy(Q_.row(j, c)); // copy current row of Q_ in the range c in R_
    r.incFirst(1);        // decrease the range
  }
}


} // namespace STK

#endif
// STK_QR_H

