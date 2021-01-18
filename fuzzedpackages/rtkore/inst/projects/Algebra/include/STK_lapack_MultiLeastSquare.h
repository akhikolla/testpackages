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
 * created on: 1 avr. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_lapack_MultiLeastSquare.h
 *  @brief In this file we define the class MultiLeastSQquare using lapack.
 **/


#ifndef STK_LAPACK_MULTILEASTSQUARE_H
#define STK_LAPACK_MULTILEASTSQUARE_H

#include "STK_ILeastSquare.h"
#include "STK_lapack_Util.h"

#include <Arrays/include/STK_CArray.h>
#include <Arrays/include/STK_CArrayVector.h>

#ifdef STKUSELAPACK

extern "C"
{

}

#endif // STKUSELAPACK


namespace STK
{

namespace lapack
{
// forward declaration
template<class ArrayB, class ArrayA> class MultiLeastSquare;
}

namespace hidden
{
/** @ingroup hidden
 *  Specialization for the lapack::MultiLeastSquare class.
 **/
template<class ArrayB_, class ArrayA_>
struct AlgebraTraits< lapack::MultiLeastSquare<ArrayB_, ArrayA_> >
{
  typedef ArrayA_ ArrayA;
  typedef ArrayB_ ArrayB;
};

} // namespace hidden


namespace lapack
{
/** @ingroup Algebra
 *  @brief The class MultiLeastSQquare solve the least square problem when
 *  the response @e b is multidimensional.
 *
 *  The class MultiLeastSquare allows to solve the least-square problem
 *  \f[
 *  \min_{x} \|b - A*x\|^2.
 *  \f]
 *  - A is an m-by-n matrix which may be rank-deficient.
 *  - B can be a vector or a Matrix.
 *
 *  It computes the minimum-norm solution to a real linear least squares
 *  problem: minimize 2-norm(| b - A*x |) using the singular value
 *  decomposition (SVD) of A.
 *  A is an M-by-N matrix which may be rank-deficient.
 *
 *  @sa STK::ILeastSquare, STK::MultiLeastSquare, STK::lapack::MultiLeastSquare
 *
 **/
template<class ArrayB, class ArrayA>
class MultiLeastSquare: public ILeastSquare< MultiLeastSquare<ArrayB, ArrayA> >
{
  public:
    typedef ILeastSquare< MultiLeastSquare<ArrayB, ArrayA> > Base;
    using Base::b_;
    using Base::a_;
    using Base::x_;
    using Base::rank_;
    /** @brief constructor
     *  @param b,a the left hand side and the right hand side of the least square problem.
     *  @param isBref,isAref are the left hand side and the right hand side references ?
     */
    MultiLeastSquare( ArrayB const& b, ArrayA const& a, bool isBref=false, bool isAref=false)
                    : Base(b, a, isBref, isAref), rcond_(-1) {};
    /** @brief template constructor
     *  @param b,a the left hand side and the right hand side of the least square
     *  problem.
     */
    template<class ArrayB_, class ArrayA_>
    MultiLeastSquare( ExprBase<ArrayB_> const& b, ExprBase<ArrayA_> const& a)
                    : Base(b, a), rcond_(-1) {}
    /** Destructor */
    virtual ~MultiLeastSquare() {};
    /** @return the condition number */
    inline Real rcond() const { return rcond_;}
    /** return the array with the singular values of A */
    inline CVectorX const& s() const { return s_;}
    /** @param rcond the condition number. If rcond<0, the machine precision is
     *  used.*/
    inline void setRcond(Real rcond) { rcond_ =rcond;}
    /** @brief solve the multi-linear least square problem.
     *  Launch gelsd LAPACK routine to perform the  decomposition.
     *  @return @c true if no error occur, @c false otherwise
     */
    bool runImpl();
    /** @brief solve the weighted least square problem.
     *  Launch gelsd LAPACK routine to perform the  decomposition.
     *  @return @c true if no error occur, @c false otherwise
     */
    template<class Weights>
    bool runImpl(Weights const& weights);

  protected:
    /** condition number used for determining the effective rank of A */
    Real rcond_;
    /** Array of the singular values */
    CVectorX s_;

  private:
    /** private method for computing the LS solution */
    bool computeLS(CArrayXX& b, CArrayXX& a);
};

/* @brief Run LS solution
 *  Launch gelsd LAPACK routine to compute the solution.
 *  @return @c true if no error occur, @c false otherwise
 */
template<class ArrayB, class ArrayA>
bool MultiLeastSquare<ArrayB, ArrayA>::runImpl()
{
  Range brows = (b_.sizeRows()< a_.sizeCols()) ? a_.cols() : b_.rows();
  // local arrays, b is resized if necessary
  CArrayXX a(a_), b(brows, b_.cols());
  b.sub(b_.rows(), b_.cols()) = b_;
  // start
  return computeLS(b,a);
}
/* @brief Run LS solution
 *  Launch gelsd LAPACK routine to compute the solution.
 *  @return @c true if no error occur, @c false otherwise
 */
template<class ArrayB, class ArrayA>
template<class Weights>
bool MultiLeastSquare<ArrayB, ArrayA>::runImpl(Weights const& w)
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
  Range brows = (b_.sizeRows()< a_.sizeCols()) ? a_.cols() : b_.rows();
  // local arrays, b is resized if necessary
  CArrayXX a(w.sqrt().diagonalize() * a_), b(brows, b_.cols());
  b.sub(b_.rows(), b_.cols()) = w.sqrt().diagonalize() * b_;
  // start
  return computeLS(b,a);
}

/* @brief Run LS solution
 *  Launch gelsd LAPACK routine to compute the solution.
 *  @return @c true if no error occur, @c false otherwise
 */
template<>
inline bool MultiLeastSquare<CArrayXX, CArrayXX>::runImpl()
{ return computeLS(b_,a_);}


/* private method for computing the LS solution */
template<class ArrayB, class ArrayA>
bool MultiLeastSquare<ArrayB, ArrayA>::computeLS(CArrayXX& b, CArrayXX& a)
{
  Range arows = a.rows(), acols = a.cols(), brows = b.rows(), bcols = b.cols();
  int m = arows.size(), n= acols.size(), nrhs = bcols.size();
  // resize if necessary
  if (brows.size()<n)
  {
    if (!b.isRef())
    {
      CArrayXX tmp(b);
      b.resize(acols,bcols);
      b.sub(brows, bcols) = tmp;
    }
    else
    { STKRUNTIME_ERROR_NO_ARG(MultiLeastSquare::computeLS,b has not enough rows);}
  }
  int lda = a.sizeRows(), ldb = b.sizeRows();
  // shift to 0
  a.shift(0,0);
  b.shift(0,0);
  s_.resize(m); s_ = 0.;
  Range srange = s_.range();
  s_.shift(0);
  // get sizes
  Real wkopt;
  int lWork = -1, iwkopt;
  int info = gelsd( m, n, nrhs
                  , a.p_data(), lda
                  , b.p_data(), ldb
                  , s_.p_data()
                  , &rcond_, &rank_
                  , &wkopt, lWork, &iwkopt);
  if( info != 0 )
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(lapack::MultiLeastSquare::computeLS,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(lapack::MultiLeastSquare::computeLS get,-info,error parameter);
    return false;
  }
  // get working sizes
  lWork = (int)wkopt;
  Real* work = new Real[lWork];
  int* iwork = new int[iwkopt];
  // Solve the least square problem
  info = gelsd( m, n, nrhs
              , a.p_data(), lda
              , b.p_data(), ldb
              , s_.p_data()
              , &rcond_, &rank_
              , work, lWork, iwork);
  delete[] iwork;
  delete[] work;
  if( info != 0 )
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(lapack::MultiLeastSquare::computeLS,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(lapack::MultiLeastSquare::computeLS get,-info,error parameter);
    return false;
  }
  // shift back
  a.shift(arows.begin(), acols.begin());
  b.shift(brows.begin(), bcols.begin());
  s_.shift(srange.begin());
  x_ = b.sub(acols,bcols);
  return true;
}

} // namespace lapack

} // namespace STK

#endif /* STK_LAPACK_MULTILEASTSQUARE_H */
