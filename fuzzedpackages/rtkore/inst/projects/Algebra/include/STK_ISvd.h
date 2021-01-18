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
 * Purpose:  Define The Interface ISvd Class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_ISvd.h
 *  @brief In this file we define the interface class ISvd.
 **/

#ifndef STK_ISVD_H
#define STK_ISVD_H

#include <Sdk/include/STK_IRunner.h>
#include <Sdk/include/STK_IRecursiveTemplate.h>
#include <STKernel/include/STK_Real.h>

#include "STK_Algebra_Util.h"

namespace STK
{
/** @ingroup Algebra
 *  @brief Compute the Singular Value Decomposition of an array.
 *
 *  The method take as:
 *  - input: A matrix A(nrow,ncol)
 *  - output:
 *    -# U Array (nrow,ncolU).
 *    -# D Vector (ncol)
 *    -# V Array (ncol,ncol).
 *  and perform the decomposition:
 *  - A = UDV'
 *  U can have more columns than A,
 *  and it is possible to compute some (all) vectors of Ker(A).
 **/
template<class Derived>
class ISvd: public IRunnerBase, public IRecursiveTemplate<Derived>
{
  protected:
    typedef typename hidden::AlgebraTraits<Derived>::ArrayU ArrayU;
    typedef typename hidden::AlgebraTraits<Derived>::ArrayD ArrayD;
    typedef typename hidden::AlgebraTraits<Derived>::ArrayV ArrayV;
    typedef typename ArrayU::Type Type;

    /** Default constructor
     *  @param A the matrix to decompose.
     *  @param ref if true, U_ is a reference of A.
     *  @param withU if @c true save the left Householder transforms in @c U_.
     *  @param withV if @c true save the right Householder transforms in @c V_.
     **/
    inline ISvd( ArrayU const& A, bool ref, bool withU = true, bool withV = true)
               : U_(A, ref), V_(), D_()
               , withU_(withU), withV_(withV), norm_(0), rank_(0), trace_(0), det_(0)
    {}
    /** constructor with other kind of array/expression
     *  @param A the matrix/expression to decompose.
     *  @param withU if @c true save the left Householder transforms in @c U_.
     *  @param withV if @c true save the right Householder transforms in @c V_.
     */
    template<class OtherDerived>
    inline ISvd( ArrayBase<OtherDerived> const& A, bool withU = true, bool withV = true)
               : U_(A), V_(), D_()
               , withU_(withU), withV_(withV), norm_(0), rank_(0), trace_(0), det_(0)
    {}
    /** Copy Constructor
     *  @param S the Svd to copy
     **/
    inline ISvd( ISvd const& S)
               : U_(S.U_, S.U_.isRef()), V_(S.V_), D_(S.D_)
               , withU_(S.withU_), withV_(S.withV_)
               , norm_(S.norm_), rank_(S.rank_), trace_(S.trace_), det_(S.det_)
    {}
    /** destructor. */
    inline virtual ~ISvd() {}
    /** Operator = : overwrite the ISvd with S.
     *  @param S the Svd to copy
     **/
    ISvd& operator=(const ISvd &S)
    {
      U_ = S.U_; V_ = S.V_;  D_ = S.D_;
      withU_ =  S.withU_; withV_ = S.withV_;
      norm_ = S.norm_; rank_ = S.rank_; trace_ = S.trace_; det_ = S.det_;
      return *this;
    }
    /** Finalize any operations that have to be done after the computation
     *  of the decomposition
     **/
    virtual void finalize()
    {
      // compute sign and trace
      trace_ = 0.;
      int s = 1;
      for (int i=D_.begin(); i< D_.end(); ++i )
      {
        Type value = D_[i];
        s *= sign(value);
        trace_ += value;
      }

      // compute norm_ (2-norm) and determinant
      norm_ = D_.maxElt();
      det_  = 0;
      if (D_.abs().minElt() > 0)
      { det_ = s * std::exp(D_.abs().log().sum());}

      // compute norm_ rank_
      rank_ = D_.size(); // full rank
      Type tol = norm_ * Arithmetic<Type>::epsilon();
      for (int i=D_.begin(); i< D_.end(); ++i )
      { if (std::abs(D_[i]) < tol) { rank_--;}}

    }

  public:
    // getters
    /** @return the product of the singular values*/
    inline Type det() const { return det_;}
    /** @return the sum of the singular values*/
    inline Type trace()  const { return trace_;}
    /** @return the norm of the matrix */
    inline Type norm()  const { return norm_;}
    /** @return the rank of the matrix */
    inline int rank()  const { return rank_;}
    /** @return U */
    inline ArrayU const& U() const { return U_;}
    /** @return  V */
    inline ArrayV const& V() const { return V_;}
    /** @return D */
    inline ArrayD const&  D() const { return D_;}

    /** implement the run method */
    virtual bool run()
    {
      if (U_.empty()) { return true;}
      // compute Svd decomposition
      if (!this->asDerived().runImpl()) return false;
      finalize();
      return true;
    }
    /** Set a new data set to ISvd class
     *  @note Take care that if U_ was previously a reference, it cannot be modified.
     *  @warning A^{-1}A give identity matrix if m<=n, and AA^{-1} give identity matrix
     *  if n<=m
     *
     *  @param A is the matrix to decompose.
     *  @param withU if @c true, we save the left Housolder transforms in U_.
     *  @param withV if @c true, we save the right Housolder transforms in V_.
     **/
    template<class OtherArray>
    void setData( OtherArray const& A, bool withU = true, bool withV = true)
    {
      U_ = A;           // Copy A in U_
      withU_ = withU;   // copy withU_ value
      withV_ = withV;   // copy withV_ value
      V_.resize(0,0), D_.resize(0);
    }
    /** Compute the generalized inverse of the matrix and put the result in res.
     *  @param res array with the result
     *  @return the result
     */
    template<class OtherArray>
    OtherArray& ginv(OtherArray& res) const;

  protected:
    /// U_ matrix
    ArrayU U_;
    /// V_ matrix
    ArrayV V_;
    /// Diagonal array of the singular values
    ArrayD D_;
    /// Compute U_ ?
    bool withU_;
    /// Compute V_ ?
    bool withV_;
    /** trace norm */
    Type norm_;
    /** rank */
    int rank_;
    /** trace norm */
    Type trace_;
    /** determinant */
    Type det_;

    /// @return the number of rows of U_
    inline int nrowU() const { return U_.sizeRows();}
    /// @return the number of columns of U_
    inline int ncolU() const { return U_.sizeCols();}
    /// @return the number of columns of D_
    inline int nrowD() const { return D_.sizeRows();}
    /// @return the number of columns of D_
    inline int ncolD() const { return D_.sizeCols();}
    /// @return the number of rows of V_
    inline int nrowV() const { return V_.sizeRows();}
    /// @return the number of columns of V_
    inline int ncolV() const { return V_.sizeCols();}
};

/* Compute the generalized inverse of the matrix and put the result in res.
 *  @param res array with the result
 *  @return the result
 */
template<class Derived>
template<class OtherArray>
OtherArray& ISvd<Derived>::ginv(OtherArray& res) const
{
  Type tol = Arithmetic<Type>::epsilon() * norm_;
  if(tol==0) { tol = Arithmetic<Type>::epsilon();}
  // compute D
  if (U_.cols() != D_.range())
  {
    Array2DDiagonal<Real> diag(U_.cols(), 0.0);
    for (int i= D_.begin(); i< D_.end(); ++i) { diag[i]=D_[i];}
    res = V_ * diag.safeInverse(tol) * U_.transpose();
  }
  else
  {
    res = V_ * D_.diagonalize().safeInverse(tol) * U_.transpose();
  }
  // compute and return UD^{-1}V'
  return (res);
}

} // namespace STK

#endif /* STK_ISVD_H */
