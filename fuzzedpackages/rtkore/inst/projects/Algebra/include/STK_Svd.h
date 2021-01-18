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
 * Purpose:  Define The Svd Class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Svd.h
 *  @brief In this file we define the Svd Class.
 **/

#ifndef STK_SVD_H
#define STK_SVD_H

#include <Arrays/include/STK_Array2DVector.h>
#include <Arrays/include/STK_Array2DDiagonal.h>
#include <Arrays/include/STK_Array2DSquare.h>
#include <Arrays/include/STK_Array2D_Functors.h>

#include "STK_ISvd.h"
#include "STK_Householder.h"
#include "STK_Givens.h"


#define MAX_ITER 30

namespace STK
{
// forward declaration
template<class Array> class Svd;

namespace hidden
{
/** @ingroup hidden
 *  Specialization for the Svd class.
 **/
template<class Array_>
struct AlgebraTraits< Svd<Array_> >
{
  typedef Array_ ArrayU;
  typedef ArrayDiagonalX ArrayD;
  typedef ArraySquareX ArrayV;
};

} // namespace hidden

/** Computing the bidiagonalization of M.
 *  The diagonal and the subdiagonal are stored in D and F
 *  @param M the matrix to bi-diagonalize, the matrix is overwritten
 *  with the left and right Householder vectors.
 *  The method return the estimate of the inf norm of M.
 *  @param D the element of the diagonal
 *  @param F the element of the surdiagnal
 **/
template<class Array>
static Real bidiag(const Array& M, ArrayDiagonalX& D, VectorX& F);
/** left eliminate the element on the subdiagonal of the row @c nrow
 *  @param D the diagonal of the matrix
 *  @param F the subdiagonal of the matrix
 *  @param nrow the number of the row were we want to rightEliminate
 *  @param U a left orthogonal Array
 *  @param withU true if we want to update U
 *  @param tol the tolerance to use
 **/
template<class Array>
static void leftEliminate( ArrayDiagonalX& D
                         , VectorX& F
                         , int const& nrow
                         , Array& U
                         , bool withU = true
                         , Real const& tol = Arithmetic<Real>::epsilon()
                         );
/** @ingroup Algebra
 *  @brief The class Svd compute the Singular Value Decomposition
 *  of a Array with the Golub-Reinsch Algorithm.
 *
 *  The method take as:
 *  - input: a matrix A(nrow,ncol)
 *  - output:
 *    -# U Array (nrow,ncol).
 *    -# D diagonal matrix (min(norw,ncol))
 *    -# V Array (ncol,ncol).
 *  and perform the decomposition:
 *  - A = UDV' (transpose V).
 *  U can have more columns than A,
 *  and it is possible to compute some (all) vectors of Ker(A).
 *
 *  @sa SK::ISvd, STK::lapack::Svd
 *
 **/
template<class Array>
class Svd: public ISvd<Svd<Array> >
{
  public :
    typedef ISvd< Svd<Array> > Base;
    typedef typename hidden::Traits<Array>::Col ColVector;
    typedef typename hidden::Traits<Array>::Row RowVector;
    using Base::U_;
    using Base::D_;
    using Base::V_;
    using Base::withU_;
    using Base::withV_;
    using Base::nrowU;
    using Base::ncolU;
    using Base::ncolV;
    using Base::norm_;
    using Base::rank_;

    /** Default constructor
     *  @param A the matrix to decompose.
     *  @param ref if true, U_ is a reference of A.
     *  @param withU if true, we save the left housolder transforms in U_.
     *  @param withV if true, we save the right housolder transforms in V_.
     **/
    inline Svd( Array const& A, bool ref= false, bool withU= true, bool withV= true)
              : Base(A, ref, withU, withV)
    {}
    /** constructor with other kind of array/expression
     *  @param A the matrix/expression to decompose.
     *  @param withU if @c true save the left housolder transforms in @c U_.
     *  @param withV if @c true save the right housolder transforms in @c V_.
     */
    template<class OtherArray>
    inline Svd( ArrayBase<OtherArray> const& A, bool withU = true, bool withV = true)
              : Base(A, withU, withV) {}
    /** Copy Constructor
     *  @param S the Svd to copy
     **/
    Svd( const Svd &S);
    /** destructor. */
    inline virtual ~Svd() {}
    /** Operator = : overwrite the Svd with S.
     *  @param S the Svd to copy
     **/
    Svd& operator=(const Svd &S);
    /** run the Svd */
    bool runImpl();
    /** Computing the diagonalization of a bi-diagonal matrix
     *  @param D the diagonal of the matrix
     *  @param F the subdiagonal of the matrix
     *  @param U a left orthogonal Array
     *  @param withU true if we want to update U
     *  @param V a right orthogonal Square Array
     *  @param withV true if we want to update V
     *  @param tol the tolerance to use
     **/
    static bool diag( ArrayDiagonalX& D
                    , VectorX& F
                    , Array& U
                    , ArraySquareX& V
                    , bool withU = true
                    , bool withV = true
                    , Real const& tol = Arithmetic<Real>::epsilon()
                    );
    /** right eliminate the element on the subdiagonal of the row nrow
     *  @param D the diagonal of the matrix
     *  @param F the subdiagonal of the matrix
     *  @param nrow the number of the row were we want to rightEliminate
     *  @param V a right orthogonal Square Array
     *  @param withV true if we want to update V
     *  @param tol the tolerance to use
     **/
    static void rightEliminate( ArrayDiagonalX& D
                              , VectorX& F
                              , int const& nrow
                              , ArraySquareX& V
                              , bool withV = true
                              , Real const& tol = Arithmetic<Real>::epsilon()
                              );
  private:
    /// Values of the Sub-diagonal
    VectorX F_;
    /// Svd main steps
    bool computeSvd();
    /// Compute U (if withU_ is true)
    void compU();
    /// Compute V (if withV_ is true)
    void compV();
};


/* Copy constructor */
template<class Array>
Svd<Array>::Svd( const Svd &S): Base(S), F_(S.F_) {}

/* Operator = : overwrite the Svd with S.*/
template<class Array>
Svd<Array>& Svd<Array>::operator=(const Svd &S)
{
  U_ = S.U_;
  V_ = S.V_;
  D_ = S.D_;
  withU_ = S.withU_;
  withV_ = S.withV_;
  norm_ = S.norm_;
  rank_ = S.rank_;
  return *this;
}

/* run the Svd */
template<class Array>
bool Svd<Array>::runImpl()
{
  if (!computeSvd()) return false;
  return true;
}


/* Main method for the svd computation. */
template<class Array>
bool Svd<Array>::computeSvd()
{
  // if the container is empty, there is nothing to do
  if (U_.empty())
  { rank_ = 0;
    norm_ = 0.0;
    return true;
  }
  int beginRow = U_.beginRows(), beginCol = U_.beginCols();
  // if U_ is just a copy of A, translate begin to 1
  // if U_ is a ref on A, this can generate an error
  U_.shift(1,1);
  // Bidiagonalize (U_)
  norm_ = bidiag(U_, D_, F_);
  // right householder vectors are in upper part of U_
  // We need to create V_ before rightEliminate
  if (withV_) { compV();}
  // rightEliminate last element of F_ if any
  if (nrowU() < ncolV())
  { rightEliminate(D_, F_, nrowU(), V_, withV_, norm_);}
  // If (U_) is not needed, we can destroy the storage
  if (withU_) { compU();}
  // Diagonalize
  bool error = diag(D_, F_, U_, V_, withU_, withV_, norm_);
  // The sub diagonal is now zero
  F_.resize(0,0);
  U_.shift(beginRow, beginCol);
  D_.shift(beginCol);
  V_.shift(beginCol);
  return error;
}


/* Bidiagonalization of the matrix M. */
template<class Array>
Real bidiag(const Array& M, ArrayDiagonalX& D, VectorX& F)
{
  typedef typename hidden::Traits<Array>::Col ColVector;
  typedef typename hidden::Traits<Array>::Row RowVector;
  // norm of the matrix M
  Real norm  = 0.0;
  // compute the number of iteration
  int begin_iter = M.beginCols();
  int last_iter  = M.beginCols() + std::min(M.sizeCols(), M.sizeRows()) -1;
  // Diagonal values
  D.resize(Range(begin_iter, last_iter, 0));
  // Upper diagonal values
  F.resize(Range(begin_iter-1, last_iter, 0));
  F.front() = 0.0;
  // Bidiagonalization of M
  // loop on the cols and rows
  Range rowRange0(M.rows())
      , rowRange1(Range(M.beginRows()+1, M.lastIdxRows(), 0))
      , colRange1(Range(M.beginCols()+1, M.lastIdxCols(), 0));
  for ( int iter=begin_iter ; iter<=last_iter
      ; iter++
      , rowRange0.incFirst(1)
      , rowRange1.incFirst(1)
      , colRange1.incFirst(1)
      )
  {
    // reference on the current column iter
    ColVector X( M, rowRange0, iter);
    // Left Householder
    D[iter] = house(X);
    // apply Householder to next cols
    applyLeftHouseholderVector(M.col(colRange1), X);
    // check if there is a row
    if ((iter < last_iter)||(M.sizeCols()>M.sizeRows()))
    {
      // ref on the current row iter
      RowVector P(M, colRange1, iter);
      // Right Householder
      F[iter] = house(P);
      // apply Householder to next rows
      applyRightHouseholderVector(M.row(rowRange1), P);
    }
    else
      F[iter] = 0.0;
    // Estimation of the norm of M
    norm = std::max(std::abs(D[iter])+std::abs(F[iter]), norm);
  }
  // return estimated norm
  return norm;
}


/* computation of V_ */
template<class Array>
void Svd<Array>::compV()
{
  // Construction of V_
  V_.resize(U_.cols());
  // Number of right Householder rotations
  int  niter = (ncolV()>nrowU()) ? (nrowU()) : (ncolV()-1);
  // initialization of the remaining rows and cols of V_ to Identity
  for (int iter=niter+2; iter<=ncolV(); iter++)
  {
    VectorX W(V_, V_.cols(), iter);
    W       = 0.0;
    W[iter] = 1.0;
  }

  Range range1(niter+1, ncolV(), 0), range2(niter+2, ncolV(), 0);
  for ( int iter0=niter, iter1=niter+1, iter2=niter+2; iter0>=1
      ; iter0--, iter1--, iter2--
      , range1.decFirst(1), range2.decFirst(1)
      )
  {
    // get beta and test
    Real beta = U_(iter0, iter1);
    if (beta)
    {
      // ref on the row iter1 of V_
      PointX  Vrow1(V_, range1, iter1);
      // diagonal element
      Vrow1[iter1] = 1.0+beta;
      // ref on the column iter1
      VectorX Vcol1(V_, range2, iter1);
      // get the Householder vector
      Vcol1 = RowVector(U_, range2, iter0);
      // Apply housholder to next cols
      for (int j=iter2; j<=ncolV(); j++)
      {
        Real aux;
        // ref on the column j
        VectorX Vcolj( V_, range2, j);
        // update column j
        Vrow1[j] = (aux = dot(Vcol1, Vcolj) * beta);
        for (int i= iter2; i <= ncolV(); ++i)
        { Vcolj[i]   += Vcol1[i] * aux;}
      }
      // compute the Householder vector
      Vcol1 *= beta;
    }
    else // nothing to do
    {
      V_(range2, iter1) = 0.0;
      V_(iter1, iter1)  = 1.0;
      V_(iter1, range2) = 0.0;
    }
  }
  // First column and rows
  V_(1,1) =1.0;
  V_(Range(2,ncolV(), 0),1) =0.0;
  V_(1,Range(2,ncolV(), 0)) =0.0;
}


/* computation of U_ */
template<class Array>
void Svd<Array>::compU()
{
  int niter = D_.size();            // Number of iterations
  int ncol  = std::min(nrowU(), ncolU()); // number of non zero cols of U_

  // initialization of the remaining cols of U_ to 0.0
  // put 0 to unused cols
  U_.col(Range(ncol+1, ncolU(), 0)) = 0.0;
  // Computation of U_
  for (int iter=niter, iter1=niter+1; iter>=1; iter--, iter1--)
  {
    // ref of the column iter
    ColVector X(U_, Range(iter1,nrowU(), 0), iter);
    // ref of the row iter
    RowVector P(U_, Range(iter,ncolU(), 0), iter);
    // Get beta and test
    Real beta = P[iter];
    if (beta)
    {
      // update the column iter
      P[iter] = 1.0 + beta;
      // Updating the cols iter+1 to ncolU_
      for (int j=iter1; j<=niter; j++)
      { // product of U_iter by U_j
        Real aux;
        ColVector Y(U_, Range(iter1, nrowU(), 0), j); // ref on the column j
        // U_j = aux = beta * X'Y
        P[j] = (aux = dot( X, Y) *beta);
        // U^j += aux * U^iter
        for (int i= iter1; i <= nrowU(); ++i)
        {
          Y[i] += X[i] * aux;
        }
      }
      // compute the vector v
      X *= beta;
    }
    else // U^iter = identity
    {
      P[iter] = 1.0;
      X = 0.0;
    }
    // update the column iter
    U_.col(Range(1,iter-1, 0), iter) = 0.0;
  }
}


/* eliminate the element of the subdiagonal with right rotations      */
template<class Array>
void Svd<Array>::rightEliminate( ArrayDiagonalX& D, VectorX& F, int const& nrow
                        , ArraySquareX& V, bool withV, Real const& tol
                        )
{
  // the element to eliminate
  Real z = F[nrow];
  // if the element is not 0.0
  if (std::abs(z)+tol != tol)
  {
    // column of the element to eliminate
    int ncol1 = nrow+1;
    // begin the Givens rotations
    for (int k=nrow, k1=nrow-1; k>=1 ; k--, k1--)
    {
      // compute and apply Givens rotation to the rows (k, k+1)
      Real aux, sinus, cosinus;
      Real y = D[k];
      D[k]   = (aux     = norm(y,z));
      z      = (sinus   = -z/aux) * F[k1];
      F[k1] *= (cosinus =  y/aux);
      // Update V_
      if (withV)
        rightGivens(V, ncol1, k, cosinus, sinus);
      // if 0.0 we can break now
      if (std::abs(z)+tol == tol) break;
    }
  }
  // the element is now 0
  F[nrow] = 0.0;        // is 0.0
}


/* eliminate the element of the subdiagonal with left rotations */
template<class Array>
void leftEliminate( ArrayDiagonalX& D, VectorX& F
                       , int const& nrow
                       , Array& U
                       , bool withU
                       , Real const& tol
                       )
{
  //the element to eliminate
  Real z = F[nrow];
  // if the element is not 0.0
  if (std::abs(z)+tol != tol)
  {
    // begin the Givens rotations
    for (int k=nrow+1; k <D.end(); k++)
    {
      // compute and apply Givens rotation to the rows (nrow, k)
      Real y = D[k];
      Real aux, cosinus, sinus;
      D[k]  = (aux     = norm(y,z));
      z     = (sinus   = -z/aux) * F[k];
      F[k] *= (cosinus = y/aux);
      // Update U_
      if (withU)
        rightGivens(U, nrow, k, cosinus, sinus);
      if (std::abs(z)+tol == tol) break;
    }
  }
  F[nrow] = 0.0;
}


/*  diagonalization of the bidiag matrix                              */
template<class Array>
bool Svd<Array>::diag( ArrayDiagonalX& D
              , VectorX& F
              , Array& U
              , ArraySquareX& V
              , bool withU
              , bool withV
              , Real const& tol
              )
{
  // result of the diag process
  bool error = false;
  // Diagonalization of A : Reduction of la matrice bidiagonale
  for (int end=D.lastIdx(); end>=D.begin(); --end)
  { // 30 iter max
    int iter;
    for (iter=1; iter<=MAX_ITER; iter++)
    { // if  the last element of the subdiagonale is 0.0
      // stop the iterations
      int beg;
      if (std::abs(F[end-1])+tol == tol)  { F[end-1] = 0.0; break;}
      // now F[end-1] !=0
      // if  D[end] == 0, we can annulate F[end-1]
      // with rotations of the columns.
      if (std::abs(D[end])+tol == tol)
      {
        D[end]   = 0.0;
        rightEliminate(D, F, end-1, V, withV, tol);
        break; // Last element of the subdiagonal is 0
      }
      // now D[end] != 0 and F[end-1] != 0
      // look for the greatest matrix such that all the elements
      // of the diagonal and subdiagonal are not zeros
      for (beg = end-1; beg>D.begin(); --beg)
      {
        if ((std::abs(D[beg])+tol == tol)||(std::abs(F[beg])+tol == tol))
        break;
      }
      // now F[beg-1]!=0
      // if D[beg] == 0 and F[beg] != 0,
      // we can eliminate the element F[beg]
      // with rotations of the rows
      if ((std::abs(D[beg])+tol == tol) && (std::abs(F[beg])+tol != tol))
      {
        D[beg] = 0.0;
        leftEliminate(D, F, beg, U, withU, tol);
      }

      // Si F[beg]==0, on augmente beg
      if (std::abs(F[beg])+tol == tol) { F[beg] = 0.0; beg++;}

      // On peut commencer les rotations QR entre les lignes beg et end
      // Shift computation
      // easy shift : commented
      // Real aux = norm(D[end],F[end-1]);
      // Real y   = (D[beg]+aux)*(D[beg]-aux);
      // Real z   = D[beg]*F[beg];
      // Wilkinson shift : look at the doc
      Real dd1 = D[end-1];      // d_1
      Real dd2 = D[end];        // d_2
      Real ff1 = F[end-2];      // f_1
      Real ff2 = F[end-1];      // f_2
      // g
      Real aux = ( (dd1+dd2)*(dd1-dd2)
                 + (ff1-ff2)*(ff1+ff2))/(2*dd1*ff2);
      // A - mu
      Real d1 = D[beg];
      Real f1 = F[beg];
      Real y  = (d1+dd2)*(d1-dd2)
              + ff2*(dd1/(aux+sign(aux,norm(1.0,aux)))- ff2);
      Real z  = d1*f1;
      // chase no null element
      int k, k1;
      for (k=beg, k1 = beg+1; k<end; ++k, ++k1)
      { // Rotation colonnes (k,k+1)
        // Input :  d1 contient D[k],
        //          f1 contient F[k],
        //          y  contient F[k-1]
        // Output : d1 contient F[k],
        //          d2 contient D[k+1],
        //          y  contient D[k]
        Real cosinus=1., sinus=0.;
        Real d2 = D[k1];
        F[k-1] = (aux = norm(y,z));                            // F[k-1]
       // arbitrary rotation if y = z = 0.0
        if (aux)
          y  = (cosinus = y/aux) * d1 - (sinus = -z/aux) * f1; // D[k]
        else
          y  = cosinus * d1 - sinus * f1;                      // D[k]

        z      = -sinus * d2;                                  // z
        d1     =  sinus * d1 + cosinus * f1;                   // F[k]
        d2    *=  cosinus;                                     // D[k+1]
        // Update V
        if (withV)
          rightGivens(V, k1, k, cosinus, sinus);
        // avoid underflow
        // Rotation lignes (k,k+1)
        // Input :  d1 contient F[k],
        //          d2 contient D[k+1],
        //          y  contient D[k]
        // Output : d1 contient D[k+1],
        //          f1 contient F[k+1],
        //          y  contient F[k]
        f1   = F[k1];
        D[k] = (aux = norm(y,z));                              // D[k]
        // arbitrary rotation if y = z = 0.0
        if (aux)
          y  = (cosinus = y/aux) * d1 - (sinus = -z/aux) * d2; // F[k]
        else
          y   = cosinus * d1 - sinus * d2;                     // F[k]

        z   = -sinus * f1;                                    // z
        d1  = sinus *d1 + cosinus * d2;                       // D[k+1]
        f1 *= cosinus;                                        // F[k+1]
        // Update U
        if (withU)
          rightGivens(U, k1, k, cosinus, sinus);
      } // end of the QR updating iteration
      D[end]   = d1;
      F[end-1] = y;
      F[beg-1] = 0.0;  // F[beg-1] is overwritten, we have to set 0.0
    } // iter
    // too many iterations
    if (iter >= 30) { error = true;}
    // positive singular value only
    if (D[end]< 0.0)
    {
      // change sign of the singular value
      D[end]= -D[end];
      // change sign of the column end of V
      if (withV)
      {
        for (int i= V.beginRows(); i <= V.lastIdxRows(); ++i)
        { V(i,end) = -V(i,end);}
      }
    }

    // We have to sort the singular value : we use a basic strategy
    Real z = D[end];        // current value
    for (int i=end+1; i<=D.lastIdx(); i++)
    { if (D[i]> z)                // if the ith singular value is greater
      { D.swap(i-1, i);       // swap the cols
        if (withU) U.swapCols(i-1, i);
        if (withV) V.swapCols(i-1, i);
      }
      else break;
    } // end sort
  } // boucle end
  return error;
}

} // namespace STK

#undef MAX_ITER

#endif
// STK_SVD_H
