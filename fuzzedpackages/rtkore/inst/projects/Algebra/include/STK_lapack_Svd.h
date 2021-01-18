/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR a PARTICULAR PURPOSE.  See the
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

/** @file STK_lapack_Svd.h
 *  @brief In this file we define the enclosing class of the dgeqrf lapack routine.
 **/


#ifndef STK_LAPACK_SVD_H
#define STK_LAPACK_SVD_H

#include <Arrays/include/STK_CArray.h>
#include <Arrays/include/STK_CArrayVector.h>

#include "STK_ISvd.h"
#include "STK_lapack_Util.h"
#include "STK_transpose.h"

namespace STK
{

namespace lapack
{
class Svd;
}

namespace hidden
{
/** @ingroup hidden
 *  Specialization for the Svd class.
 **/
template<>
struct AlgebraTraits< lapack::Svd >
{
  typedef CArrayXX ArrayU;
  typedef CVectorX ArrayD;
  typedef CArrayXX ArrayV;
};

} // namespace hidden


namespace lapack
{
/** @ingroup Algebra
 *  {
 *    @class Svd
 *    @brief Svd computes the SVD decomposition of a real matrix using the
 *    Lapack routine dgeqrf.
 *
 *  The method take as:
 *  - input: a matrix A(nrow,ncol)
 *  - output:
 *    -# U Array (nrow,ncol).
 *    -# D diagonal matrix (min(norw,ncol))
 *    -# V Array (ncol,ncol).
 *  and perform the decomposition:
 *  - A = UDV' (transpose V).
 *
 *    Output can be tuned using @c jobz option:
 *    @verbatim
 *      jobz is Char*1
 *      Specifies options for computing all or part of the matrix u:
 *      = 'A':  all m columns of u and all n rows of v**T are
 *              returned in the arrays u and vt;
 *      = 'S':  the first min(m,n) columns of u and the first
 *              min(m,n) rows of v**T are returned in the arrays U
 *              and vt;
 *      = 'O':  If m >= n, the first n columns of u are overwritten
 *              on the array a and all rows of v**T are returned in
 *              the array vt;
 *              otherwise, all columns of u are returned in the
 *              array u and the first m rows of v**T are overwritten
 *              in the array a;
 *      = 'N':  no columns of u or rows of v**T are computed.
 *   @endverbatim
 *
 * @sa SK::ISvd, STK::Svd
 **/
class Svd: public ISvd<Svd>
{
  public:
    typedef ISvd< Svd > Base;
    typedef typename hidden::Traits<CArrayXX>::Col ColVector;
    typedef typename hidden::Traits<CArrayXX>::Row RowVector;

    using Base::norm_;
    using Base::rank_;
    using Base::U_;
    using Base::D_;
    using Base::V_;
    using Base::withU_;
    using Base::withV_;
    using Base::nrowU;
    using Base::ncolU;
    using Base::ncolV;

    /** Default constructor.
     *  @param A matrix to decompose
     *  @param ref if @c true, U_ is a reference of A.
     *  @param withU if @c true save the left Householder transforms in @c U_.
     *  @param withV if @c true save the right Householder transforms in @c V_.
     **/
    inline Svd( CArrayXX const&  A, bool ref = false, bool withU= true, bool withV= true)
              : Base(A, ref, withU, withV), jobz_( (withU|withV) ? 'O':'N') {}
    /** constructor with other kind of array/expression
     *  @param A the matrix/expression to decompose.
     *  @param withU if @c true save the left Householder transforms in @c U_.
     *  @param withV if @c true save the right Householder transforms in @c V_.
     */
    template<class OtherArray>
    inline Svd( ArrayBase<OtherArray> const& A, bool withU= true, bool withV= true)
              : Base(A, withU, withV), jobz_( (withU|withV) ? 'O':'N') {}
    /** Copy constructor.
     *  @param decomp the decomposition to copy
     **/
    inline Svd( Svd const& decomp): Base(decomp), jobz_(decomp.jobz_) {}
    /** virtual destructor */
    inline virtual ~Svd() {}
    /** @brief clone pattern */
    inline virtual Svd* clone() const { return new Svd(*this);}
    /** Operator = : overwrite the Svd with decomp. */
    inline Svd& operator=(Svd const& decomp)
    {
      Base::operator=(decomp);
      jobz_ = decomp.jobz_;
      return *this;
    }
    // getters
    /** @return the option chosen for the svd */
    char jobz() const { return jobz_;}
    // setters
    /** set the option chosen for the svd */
    void setJobz(char jobz) { jobz_ = jobz;}

    /** @brief Run svd decomposition */
    bool runImpl();

  private:
     /** option */
     char jobz_;
     /** compute the svd decomposition. a contains either u, vt (if jobz_=='O')
      *  or is destroyed at the end of the oputput. */
     bool computeSvd(CArrayXX& a, CArrayXX& u, CVectorX& s, CArrayXX& v);
};

/* @brief Run svd decomposition */
inline bool Svd::runImpl()
{
  if (jobz_ == 'A' || jobz_ == 'S')
  {
    msg_error_ = _T("In lapack::Svd::runImpl, the options 'A' and 'S' are not available");
    return false;
  }
  int beginRow = U_.beginRows(), beginCol = U_.beginCols();
  bool result = true;
  // compute results
  if ( (jobz_ == 'O' && U_.sizeRows() >= U_.sizeCols()) || jobz_ == 'N')
  {
    if (!computeSvd(U_, U_, D_, V_)) { result = false;}
  }
  else
  { // jobz_ == 'O' and m<n, V_ will contain u and U_ will contain vt
    if (!computeSvd(U_, V_, D_, V_)) { result = false;};
    U_.exchange(V_); // !now U_ is (m,m) and V_ is (m,n)
  }
  U_.shift(beginRow, beginCol);
  D_.shift(beginCol);
  transpose(V_.shift(beginCol));
  return result;
}


/* Computing the Svd decomposition of the matrix U_. */
inline bool Svd::computeSvd(CArrayXX& a, CArrayXX& u, CVectorX& s, CArrayXX& vt)
{
  int m = a.sizeRows(), n = a.sizeCols(), nbSv = std::min(m,n);
  a.shift(0,0);
  // Workspace and status variables:
  double workSize;
  double *p_work = &workSize;
  int* p_iwork = new int[8*nbSv];
  int lwork = -1;
  int info;
  //
  // Call dgesdd_ with lwork = -1 to query optimal workspace size:
  info = gesdd( jobz_, m, n
               , a.p_data(), m, s.p_data(), u.p_data(), m, vt.p_data(), n
               , p_work, lwork, p_iwork);
  // check any error
  if (info!=0)
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(lapack::Svd::computeSvd,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(lapack::Svd::computeSvd get,-info,error parameter);
    return false;
  }
  // optimal storage dimension is returned in work[0]
  lwork = workSize;
  p_work = new Real[lwork];
  // resize u
  if ( !((jobz_ == 'O' && m >= n) || (jobz_ == 'N')) )
  {
    int ucol = (jobz_ == 'A' || jobz_ == 'O') ? m : nbSv;
    u.resize(m, ucol);
  }
  if ( !((jobz_ == 'O' && m < n)  || (jobz_ == 'N')) )
  {
    int ldvt = (jobz_ == 'A' || jobz_ == 'O') ? n : nbSv;
    vt.resize(ldvt,n);
  }
  s.resize(nbSv).shift(0);
  u.shift(0,0);
  vt.shift(0,0);

  // Call dgesdd_ to do the actual computation:
  info = gesdd( jobz_, m, n
              , a.p_data(), a.sizeRows(), s.p_data()
              , u.p_data(), u.sizeRows()
              , vt.p_data(), vt.sizeRows()
              , p_work, lwork, p_iwork);
  // clean
  delete[] p_work;
  delete[] p_iwork;
  // check any error
  if (info!=0)
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(lapack::Svd::computeSvd,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(lapack::Svd::computeSvd get,-info,error parameter);
    return false;
  }
  return true;
}

} // namespace lapack

} // namespace STK

#endif /* STK_LAPACK_SVD_H */
