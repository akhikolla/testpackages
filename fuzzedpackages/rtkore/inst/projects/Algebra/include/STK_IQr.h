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
 * Purpose:  Define The Interface SymEigen Class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IQr.h
 *  @brief In this file we define the Interface class IQr (QR decomposition).
 **/

#ifndef STK_IQR_H
#define STK_IQR_H

#include <STKernel/include/STK_Real.h>
#include <Sdk/include/STK_IRunner.h>

#include <Arrays/include/STK_Array2DVector.h>
#include <Arrays/include/STK_Array2DUpperTriangular.h>

#include "STK_Algebra_Util.h"
#include "STK_Householder.h"
#include "STK_Givens.h"

#ifdef STK_ALGEBRA_DEBUG
#include "Arrays/include/STK_Display.h"

template< class Container2D >
void print(Container2D const& A, STK::String const& name)
{
  stk_cout << "print: " << name << _T("\n";);
  stk_cout << name << _T(".isRef() =")        << A.isRef()  << _T("\n");
  stk_cout << name << _T(".cols() =")      << A.cols()  << _T("\n");
  stk_cout << name << _T(".rows() =")      << A.rows()  << _T("\n\n");
  stk_cout << name << _T(".rangeCols().isRef() =")  << A.rangeCols().isRef() << _T("\n");
  stk_cout << name << _T(".rangeCols() =\n")  << A.rangeCols() << _T("\n");
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
/** @ingroup Algebra
 *  @brief The class IQr is an interface class for the method
 *  computing the QR Decomposition of a ArrayXX.
 *
 *  The QR decomposition of a matrix require
 *  - Input:  A matrix of size M-by-N
 *  - Output:
 *     -# Q Array of size M-by-N.
 *     -# R Upper Triangular matrix of dimension min(M,N)-by-N
 *     -# \f$ A = QR \f$
 *
 *  Computation of matrices @c Q and @c R is triggered by the @c run() method.
 *  Some tuning parameters can be adjusted before running.
 *
 *  The matrix Q is represented as a product of elementary reflectors
 *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
 *
 *  Each H(i) has the form
 *     H(i) = I - tau * v * v'
 *
 *  @note This interface class implement the virtual @c run method of @c IRunnerBase
 *  by calling the derived method @c runImpl
 *
 *
 **/
template <class Derived>
class IQr: public IRunnerBase, public IRecursiveTemplate<Derived>
{
  protected:
    typedef typename hidden::AlgebraTraits<Derived>::Array Array;
    typedef typename hidden::AlgebraTraits<Derived>::ColVector ColVector;
    typedef typename hidden::AlgebraTraits<Derived>::RowVector RowVector;
    /** @brief Constructor
     *  @param Q The matrix to decompose
     *  @param ref is Q_ a reference ? (avoid data copy, but the original data
     *  set is destroyed)
     */
    inline IQr( Array const& Q, bool ref = false)
              : IRunnerBase(), Q_(Q, ref), compq_(false), rank_(0)
    {
      if (Q.beginRows() != Q.beginCols())
        STKRUNTIME_ERROR_NO_ARG(IQR::IQR,Wrong data set: Q.beginRows() must be equal to Q.beginCols());
    }
    /** @brief Constructor
     *  @param Q The matrix to decompose
     */
    template<class DerivedExpr>
    IQr( ExprBase<DerivedExpr> const& Q): IRunnerBase(), compq_(false), rank_(0)
    {
      if (Q.beginRows() != Q.beginCols())
      { STKRUNTIME_ERROR_NO_ARG(IQR::IQR,Wrong data set: Q.beginRows() must be equal to Q.beginCols());}
      Q_ = Q.asDerived();
    }
    /** Copy constructor.
     *  @param decomp the decomposition to copy
     **/
    inline IQr( IQr const& decomp)
              : Q_(decomp.Q_), R_(decomp.R_), compq_(decomp.compq_), rank_(decomp.rank_)
    {}

  public :
    /** virtual destructor */
    inline virtual ~IQr() {}
    /** Operator = : overwrite the this with decomp. */
    IQr& operator=(IQr const& decomp)
    {
      Q_ = decomp.Q_;
      R_ = decomp.R_;
      compq_ = decomp.compq_;
      rank_ = decomp.rank_;
      return *this;
    }
    // getters
    /** @return @c true if Q_ is computed, @c false otherwise */
    inline bool isCompQ() const { return compq_;}
    /** @return the rank of the matrix Q */
    inline int rank() const { return rank_;}
    /** give the matrix Q of the QR decomposition.
     *  @return the matrix Q of the QR decomposition
     **/
    inline Array const& Q() const  { return Q_;}
    /** give the matrix R of the QR decomposition.
     *  @return the matrix R of the QR decomposition
     **/
    inline ArrayUpperTriangularXX const& R() const { return R_;}

    // setters
    /** Set the compq_ flag. By default compq is @c true.
     *  @param compq @c true if we want to compute the @c Q matrix, @c false otherwise
     **/
    inline void setCompQ(bool compq) { compq_ = compq;}

    // methods
    /** Compute the QR decomposition.
     *  Delegate to Derived classes the concrete computation of the
     *  decomposition using @c runImpl method.
     *  @return @c true if the computation succeed, @c false otherwise
     **/
    virtual bool run();
    /** Compute Q (to use after run). After the run process, Q_ store
     *  the householder vector in its column. Call compQ, if you want to
     *  obtain Q in its true form.
     *  Without effect if (compq_ == true)
     **/
    void compQ();
    /** Delete the n last columns and update the QR decomposition.
     *  @param n number of column to delete
     **/
    void popBackCols(int n =1);
    /** Delete the column pos and update the QR decomposition.
     *  @param pos the position of the column to delete
     **/
    void eraseCol(int pos);
    /** Add a column with value T and update th QR decomposition.
     *  @param T the column to add
     **/
    template<class ColVector_>
    void pushBackCol(ColVector_ const& T);
    /** Add a column with value T at the position pos and update the QR
     *  decomposition.
     *  @param T the column to insert
     *  @param pos the position of the column to insert
     **/
    template<class ColVector_>
    void insertCol(ColVector_ const& T, int pos);
    /** overloading of setData.
     * @param data the data set to set.
     **/
    template<class OtherDerived>
    void setData( ExprBase<OtherDerived> const& data)
    { Q_ = data.asDerived(); R_.clear(); compq_ = false; rank_ = 0;}

    /* TODO : Delete the ith row and update the QR decomposition
     * TODO : Add a row with value T and update th QR decomposition
     */
    //void popBackRows();
    //void eraseRows(int pos);
    //void pushBackRows(RowVector const& T);
    //void insertRows(RowVector const& T, int pos);

  protected :
    /** Q Array of the QR decomposition */
    Array Q_;
    /** R Array of th QR decomposition */
    ArrayUpperTriangularXX R_;
    /// is Q computed ?
    bool compq_;
    /// estimated rank_ of the matrix A
    int rank_;
};

template<class Derived>
bool IQr<Derived>::run()
{
  // if Q is empty, default values are ok
  if (Q_.empty()) { compq_ = true; return true;}
  // compute QR decomposition
  if (!this->asDerived().runImpl()) return false;
  // compute Q if asked
  if (compq_) compQ();
  // compute rank_ using maximal number of R_
  int r = R_.beginRows();
  Real tol =  Arithmetic<Real>::epsilon() * std::abs(R_(r,r));

  int n = std::min(R_.endRows(), R_.endCols());
  for(r= R_.beginRows(), rank_=0; r<n; ++r, ++rank_) { if (std::abs(R_(r,r)) < tol) break;}
  return true;
}

template<class Derived>
void IQr<Derived>::compQ()
{
#ifdef STK_ALGEBRA_VERBOSE
  stk_cout << _T("Entering IQr::compQ()") << _T("\n");
#endif
  // if Q_ is computed yet
  if (compq_) return;
  // number of non zero cols of Q_
  int ncol  = std::min(Q_.sizeRows(), Q_.sizeCols()), lastCol;
  // add or remove the column
  if (ncol < Q_.sizeCols())
  {
    Q_.popBackCols(Q_.sizeCols() - ncol);
    Q_.popBackCols(Q_.sizeCols() - ncol);
    lastCol = Q_.lastIdxCols();
  }
  else
  {
    lastCol = Q_.lastIdxCols();
    if (ncol < Q_.sizeRows())
    {
      Q_.pushBackCols(Q_.sizeRows() -ncol);
      // Initialize added columns
      Q_.col( _R( lastCol+1, Q_.lastIdxCols()) ).setValue(0);
      for (int i=lastCol+1; i< Q_.endCols(); ++i) { Q_(i, i) = 1.0;}
    }
  }
  // compute other columns
  for (int i=lastCol; i>=Q_.beginCols(); i--)
  {
    // create reference on current householder vector
    ColVector u(Q_, _R(i, Q_.lastIdxRows()), i);
    // Apply Householder vector to the right of the matrix
    applyLeftHouseholderVector( Q_( _R(i, Q_.lastIdxRows()), _R(i+1, Q_.lastIdxCols())), u);
    // update the ith column
    Q_( _R(Q_.beginRows(),i-1)   , i ) = 0.0;   //     0:(i-1)
    Q_( _R(i+1, Q_.lastIdxRows()), i ) *= Q_(i,i); // (i+1):M
    Q_( i                        , i ) += 1.0;  //     i:i
    // update the column i
  }
  // Q_ is now computed
  compq_ = true;
#ifdef  STK_ALGEBRA_VERBOSE
  stk_cout << _T("Terminating IQr::compQ().") << _T("\n");
#endif
}

template<class Derived>
void IQr<Derived>::popBackCols(int n) { R_.popBackCols(n);}

/* Delete the jth column and update the QR decomposition
 **/
template<class Derived>
void IQr<Derived>::eraseCol(int pos)
{
  if (pos < R_.beginCols())
  { STKOUT_OF_RANGE_1ARG(Qr::eraseCol,pos,pos<R_.beginCols());}
  if (R_.lastIdxCols() < pos)
  { STKOUT_OF_RANGE_1ARG(Qr::eraseCol,pos,pos<R_.lastIdxCols()<pos);}
  // if Q_ is not computed yet
  if (!compq_) compQ();
  // compute the number of iteration for updating to zeroed
  int niter = std::min(R_.lastIdxCols(), R_.lastIdxRows());//R_.beginCols()-1+std::min(R_.sizeRows(), R_.sizeCols());
  // Zeroed the remaining elements (z)
  for (int iter = pos+1; iter<= niter; iter++)
  {
    Real sinus, cosinus;
    // compute the Givens rotation
    R_(iter-1, iter) = compGivens( R_(iter-1, iter), R_(iter, iter), cosinus, sinus);
    R_(iter, iter)   = 0.0;
    // if necessary update R_ and Q_
    if (sinus)
    {
      // create a reference on the sub-ArrayXX
      ArrayUpperTriangularXX Rsub(R_.col( _R(iter+1, R_.lastIdxCols()) ), true);
      // Update the next rows (iter1:ncolr_) of R_
      leftGivens(Rsub, iter-1, iter, cosinus, sinus);
      // Update the cols of Q_
      rightGivens(Q_, iter-1, iter, cosinus, sinus);
    }
  }
  // erase the column pos
  R_.eraseCols(pos);
  // update the range of the remaining cols of the container
  R_.update( Range(pos, std::min(R_.lastIdxRows(), R_.lastIdxCols()), 0));
}

/* Adding the last column and update the QR decomposition. */
template<class Derived>
template<class ColVector_>
void IQr<Derived>::pushBackCol(ColVector_ const& T)
{
  STK_STATIC_ASSERT(ColVector_::structure_==(int)Arrays::vector_||ColVector_::structure_==(int)Arrays::point_,YOU_HAVE_TO_USE_A_VECTOR_OR_POINT_IN_THIS_METHOD)
  // check conditions
  if (T.range() != Q_.rows())
  { STKRUNTIME_ERROR_NO_ARG(Qr::pushBackCol,T.range() != Q_.rows());}
  // if Q_ is not computed yet
  if (!compq_) compQ();
  // Adding a column to R
  int lastColR = R_.endCols();
  // Create an auxiliary container
  VectorX Rncolr = Q_.transpose() * T; // Rncolr of size Q_.cols()
  // update Q_
  for (int iter = Q_.lastIdxCols()-1, iter1 = Q_.lastIdxCols(); iter>=lastColR; iter--, iter1--)
  {
    Real sinus, cosinus;
    // compute the Givens rotation
    Rncolr[iter] = compGivens( Rncolr[iter], Rncolr[iter1], cosinus, sinus);
    // apply Givens rotation if necessary
    if (sinus)
    { rightGivens(Q_, iter, iter1, cosinus, sinus);}
  }
  // update R_
  R_.pushBackCols();
  R_.col(lastColR).copy(Rncolr.sub(R_.rangeRowsInCol(lastColR)));
}


/* Add a column with value T at the position pos and update the QR
 *  decomposition.
 *  @param T the column to insert
 *  @param pos the position of the column to insert
 **/
template<class Derived>
template<class ColVector_>
void IQr<Derived>::insertCol(ColVector_ const& T, int pos)
{
  STK_STATIC_ASSERT(ColVector_::structure_==(int)Arrays::vector_||ColVector_::structure_==(int)Arrays::point_,YOU_HAVE_TO_USE_A_VECTOR_OR_POINT_IN_THIS_METHOD)
  if (pos < R_.beginCols())
  { STKOUT_OF_RANGE_1ARG(Qr::insertCol,pos,pos<R_.beginCols());}
  if (R_.lastIdxCols() < pos)
  { STKOUT_OF_RANGE_1ARG(Qr::insertCol,pos,pos<R_.lastIdxCols()<pos);}
  if (T.range() != Q_.rows())
  { STKRUNTIME_ERROR_1ARG(Qr::insertCol,pos,T.range() != Q_.rows());}
  // if Q_ is not computed yet
  if (!compq_) compQ();
  // Adding a column to R
  R_.insertCols(pos);
  // update the range of the remaining cols of R_
  R_.update( _R(pos+1, std::min(R_.lastIdxRows(), R_.lastIdxCols())) );
  for (int i=pos+1; i< std::min(R_.endRows(), R_.endCols()); ++i) R_(i,i) = 0.0;

  VectorX Rpos =  Q_.transpose() * T;
  // Zeroed the unwanted elements
  for (int iter= Q_.lastIdxCols(), iterm1= Q_.lastIdxCols()-1; iter>pos; iterm1--, iter--)
  {
    Real sinus, cosinus;
    // compute the Givens rotation
    Rpos[iterm1]  = compGivens(Rpos[iterm1], Rpos[iter], cosinus, sinus);
    // apply Givens rotation if necessary
    if (sinus)
    {
      // create a reference on the sub-ArrayXX
      ArrayUpperTriangularXX Rsub(R_.col(_R(iter, R_.lastIdxCols())), true);
      // Update the next rows (iter:ncolr_) of R_
      leftGivens( Rsub, iterm1, iter, cosinus, sinus);
      // Update the cols of Q_
      rightGivens(Q_, iterm1, iter, cosinus, sinus);
    }
  }
  // update R_
  R_.col(pos) = Rpos.sub(R_.rangeRowsInCol(pos));
  R_.update(pos);
}

} // namespace STK

#endif //STK_IQR_H
