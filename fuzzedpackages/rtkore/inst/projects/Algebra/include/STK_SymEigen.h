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
 * Purpose:  Define The SymEigen Class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_SymEigen.h
 *  @brief In this file we define the SymEigen class (for a symmetric matrix).
 **/

#ifndef STK_SYMEIGEN_H
#define STK_SYMEIGEN_H

#include "STK_ISymEigen.h"
#include "Arrays/include/STK_Array2DVector.h" // for F_
#include "STK_Householder.h"
#include "STK_Givens.h"

#define MAXITER 30

namespace STK
{
// forward declaration
template<class SquareArray> class SymEigen;

namespace hidden
{
/** @ingroup hidden
 *  Specialization for the Qr class.
 **
 **/
template<class SquareArray_>
struct AlgebraTraits< SymEigen<SquareArray_> >
{
  typedef SquareArray_ SquareArray;
};

} // namespace hidden

/** @ingroup Algebra
 *  @brief The class SymEigen compute the eigenvalue Decomposition
 *  of a symmetric ArrayXX.
 *
 *  The decomposition of a symmetric matrix require
 *  - Input:  A symmetric matrix A of size (n,n)
 *  - Output:
 *     -# P Array of size (n,n).
 *     -# D Vector of dimension n
 *     -# \f$ A = PDP' \f$
 *  The matrix A can be copied or overwritten by the class.
 *
 *  The 2-norm (operator norm) of the matrix is given. if the 2-norm is less
 *  than the arithmetic precision of the type @c Real, the rank is not
 *  full.
 *  Thus the user can be faced with a deficient rank matrix and with a norm and
 *  a determinant very small (but not exactly 0.0).
 *
 *  @sa STK::ISymEigen, STK::lapack::SymEigen
 **/
template<class SquareArray>
class SymEigen: public ISymEigen<SymEigen<SquareArray> >
{
  public:
    typedef ISymEigen<SymEigen<SquareArray> > Base;
    using Base::eigenValues_;
    using Base::eigenVectors_;
    using Base::range_;
    using Base::norm_;
    using Base::rank_;
    using Base::det_;

    /** @brief Default Constructor */
    SymEigen();
    /** @brief Constructor
     *  @param data reference on a symmetric square matrix
     *  @param ref @c true if we overwrite the data set, @c false otherwise
     *  @note data can be a reference if and only if it is a CSquareX
     */
    SymEigen( SquareArray const& data, bool ref =false);
    /** constructor.
     *  @param data A reference on a symmetric expression matrix to decompose.
     **/
    template<class Derived>
    SymEigen( ExprBase<Derived> const& data);
    /** Copy constructor.
     *  @param S the EigenValue to copy
     **/
    SymEigen( SymEigen const& S);
    /** virtual destructor */
    virtual ~SymEigen() {}

    /** @brief Diagonalization of eigenVectors_
     *  @return @c true if no error occur, @c false otherwise
     * */
    bool runImpl();
    /** Operator = : overwrite the SymEigen with S.
     *  @param S SymEigen to copy
     *  @return a reference on this
     **/
    inline SymEigen& operator=(const SymEigen &S)
    {
      Base::operator=(S);
      return *this;
    }

  private:
    /** Temporary vector. Values of the subdiagonal. */
    VectorX F_;
    /** @brief compute the tri-diagonalization of eigenVectors_ */
    void tridiagonalize();
    /** @brief compute the Householder matrix and P */
    void compHouse();
    /** computing the diagonalization of eigenValues_ and F_ */
    void diagonalize();
};


/* @brief Default Constructor */
template<class SquareArray>
SymEigen<SquareArray>::SymEigen(): Base() {}
/* @brief Constructor
 *  @param data reference on a symmetric square matrix
 *  @param ref @c true if we overwrite the data set, @c false otherwise
 */
template<class SquareArray>
SymEigen<SquareArray>::SymEigen( SquareArray const& data, bool ref)
                               : Base(data)
{}
/* @brief Constructor
 *  @param data reference on a symmetric square matrix
 *  @param ref @c true if we overwrite the data set, @c false otherwise
 */
template<>
inline SymEigen<CSquareX>::SymEigen( CSquareX const& data, bool ref)
                                   : Base(data, ref)
{}

/* constructor.
 *  @param data A reference on the symmetric matrix to decompose.
 **/
template<class SquareArray>
template<class Derived>
SymEigen<SquareArray>::SymEigen( ExprBase<Derived> const& data): Base(data) {}
/* Copy constructor.
*  @param S the EigenValue to copy
**/
template<class SquareArray>
SymEigen<SquareArray>::SymEigen( SymEigen const& S): Base(S) {}

/* Main methods.  */
/* Compute diagonalization of the symmetric matrix */
template<class SquareArray>
bool SymEigen<SquareArray>::runImpl()
{
#ifdef STK_ALGEBRA_VERBOSE
    stk_cout << _T("Entering in SymEigen::run()\n");
#endif
  try
  {
    // copy data
    eigenValues_.resize(eigenVectors_.range());
    norm_ = 0.;
    rank_ = 0;
    det_ = 0.;
#ifdef STK_ALGEBRA_VERBOSE
    stk_cout << _T("calling SymEigen::tridiagonalize()\n");
#endif
    // tridiagonalize eigenVectors_
    tridiagonalize();
#ifdef STK_ALGEBRA_VERBOSE
    stk_cout << _T("calling SymEigen::compHouse()\n");
#endif
    // compute eigenVectors_
    compHouse();
#ifdef STK_ALGEBRA_VERBOSE
    stk_cout << _T("calling SymEigen::diagonalize()\n");
#endif
    // Diagonalize
    diagonalize();
#ifdef STK_ALGEBRA_VERBOSE
    stk_cout << _T("calling SymEigen::finalize()\n");
#endif
    // compute rank, norm and determinant
    this->finalizeStep();
  }
  catch (Exception const& e)
  {
    this->msg_error_ = e.error();
    return false;
  }
  return true;
}

/* tridiagonalization of the symmetric matrix eigenVectors_. Only the lower
 * part of eigenVectors_ used. eigenVectors_ is overwritten with the Householder
 * vectors. eigenValues_ contains the diagonal.
 **/
template<class SquareArray>
void SymEigen<SquareArray>::tridiagonalize()
{
  // Upper diagonal values
  F_.resize(Range(range_.begin()-1, range_.lastIdx(), 0));
  F_.front() = 0.0; F_.back() =0.0;
  // initial range of the Householder vectors
  Range range1(Range(range_.begin()+1, range_.lastIdx(), 0));
  Range range2(Range(range_.begin()+2, range_.lastIdx(), 0));
  // Bidiagonalization of eigenVectors_
  // loop on the cols and rows
  for ( int i=range_.begin(), i1=range_.begin()+1, i2=range_.begin()+2; i<range_.lastIdx()
      ; i++, i1++, i2++, range1.incFirst(1), range2.incFirst(1)
      )
  {
    // ref on the current column iter in the range iter1:last
    typename hidden::Traits<CSquareX>::Col v(eigenVectors_.col(range1, i));
    // Compute Householder vector and get subdiagonal element
    F_[i] = house(v);
    // Save diagonal element
    eigenValues_[i] = eigenVectors_(i,i);
    // get beta
    Real beta = v.front();
    if (beta)
    {
      // ref on the current column iter1 in the range iter1:last
      typename hidden::Traits<CSquareX>::Col M1(eigenVectors_.col(range1, i1));
      // aux1 will contain <v,p>
      Real aux1 = 0.0;
      // apply left and right Householder to eigenVectors_
      // compute D(iter1:range_.last_) = p and aux1 = <p,v>
      // Computation of p_iter1 = beta * eigenVectors_ v using the lower part of eigenVectors_
      // save p_iter1 in the unused part of eigenValues_ and compute p'v
      aux1 += (eigenValues_[i1] = beta*(M1[i1] + M1.sub(range2).dot(v.sub(range2)))) /* *1.0 */;
      // other cols
      for (int k=i2; k<range_.end(); k++)
      {
        // Computation of p_i = beta * eigenVectors_ v using the lower part of eigenVectors_
        // save p_i in the unusued part of eigenValues_ and compute p'v
        Real aux = M1[k] /* *1.0 */;
        for (int j=i2;  j<=k;     j++) { aux += eigenVectors_(k,j)*v[j];}
        for (int j=k+1; j<range_.end(); j++) { aux += eigenVectors_(j,k)*v[j];}
        // save p_i in the unusued part of eigenValues_ and compute p'v
        aux1 += (eigenValues_[k] = beta*aux) * v[k];
      }
      // update diagonal element M_ii+= 2 v_i * q_i = 2* q_i (i=iter1)
      // aux = q_iter1 and aux1 = beta*<p,v>/2 (we don't save aux in eigenValues_)
      Real aux = (eigenValues_[i1] + (aux1 *= (beta/2.0)));
      M1[i1] += 2.0*aux;
      // update lower part: all rows
      // compute q_i and update the lower part of eigenVectors_
      for (int k=i2; k<range_.end(); k++)
      {
        // get q_i and save it in eigenValues_i=q_i = p_i + <p,v> * beta * v_i/2
        eigenValues_[k] += aux1 * v[k];
        // Compute eigenVectors_ + u q' + q u',
        // update the row i, range_.begin() element
        M1[k] += eigenValues_[k] /* *1.0 */+ v[k] * aux;
        // update the row i: all cols under the diagonal
        for (int j=i2; j<=k; j++)
          eigenVectors_(k,j) += v[j] * eigenValues_[k] + v[k] * eigenValues_[j];
      }
    }
  }
  // Last col
  eigenValues_[range_.lastIdx()] = eigenVectors_(range_.lastIdx(),range_.lastIdx());
}

// Compute eigenVectors_ from the Householder rotations
template<class SquareArray>
void SymEigen<SquareArray>::compHouse()
{
  // compute eigenVectors_
  // iter0 is the column of the Householder vector
  // iter is the current column to compute
  for ( int iter0=range_.lastIdx()-1, iter=range_.lastIdx(), iter1=range_.lastIdx()+1
      ; iter0>=range_.begin()
      ; iter0--, iter--, iter1--)
  {
    // reference on the Householder vector
    typename hidden::Traits<CSquareX>::Col v(eigenVectors_.col(Range(iter, range_.lastIdx(), 0), iter0));
    // Get Beta
    Real beta = v[iter];
    if (beta)
    {
      // Compute Col iter -> eigenVectors_ e_{iter}= e_{iter}+ beta v
      eigenVectors_(iter, iter) = 1.0 + beta /* *1.0 */;
      for (int i=iter1; i<range_.end(); i++)
      { eigenVectors_(i, iter) = beta * v[i];}

      // Update the other Cols
      for (int i=iter1; i<range_.end(); i++)
      {
        // compute beta*<v, eigenVectors_^i>
        Real aux = 0.0;
        for (int j=iter1; j<range_.end(); j++)
        { aux += eigenVectors_(j, i) * v[j]; }
        aux *= beta;
        // range_.begin() row (iter)
        eigenVectors_(iter, i) = aux;
        // other rows
        for (int j=iter1; j<range_.end(); j++)
        { eigenVectors_(j, i) += aux * v[j];}
      }
    }
    else // beta = 0, nothing to do
    { eigenVectors_(iter,iter)=1.0;
      for (int j=iter1; j<range_.end(); j++)
      { eigenVectors_(iter, j) =0.0; eigenVectors_(j, iter) = 0.0;}
    }
  }
  // range_.begin() row and range_.begin() col
  eigenVectors_(range_.begin(), range_.begin()) = 1.0;
  for (int j=range_.begin()+1; j<range_.end(); j++)
  { eigenVectors_(range_.begin(), j) = 0.0; eigenVectors_(j, range_.begin()) = 0.0;}
}

// diagonalize eigenValues_ and F_
template<class SquareArray>
void SymEigen<SquareArray>::diagonalize()
{
  // Diagonalisation of eigenVectors_
  for (int iend=range_.lastIdx(); iend>=range_.begin()+1; iend--)
  {
    int iter;
    for (iter=0; iter<MAXITER; iter++) // fix the number of iterations max
    {
      // check cv for the last element
      Real sum = std::abs(eigenValues_[iend]) + std::abs(eigenValues_[iend-1]);
      // if the first element of the small subdiagonal
      // is 0. we stop the QR iterations and increment iend
      if (std::abs(F_[iend-1])+sum == sum)
      { F_[iend-1] = 0.0 ; break;}
      // look for a single small subdiagonal element to split the matrix
      int ibeg = iend-1;
      while (ibeg>range_.begin())
      {
        ibeg--;
        // if a subdiagonal is zero, we get a sub matrix unreduced
        //sum = std::abs(eigenValues_[ibeg])+std::abs(eigenValues_[ibeg+1]);
        if (std::abs(F_[ibeg])+std::abs(eigenValues_[ibeg]) == std::abs(eigenValues_[ibeg]))
        { F_[ibeg] = 0.; ibeg++; break;}
      };
      // QR rotations between the rows/cols ibeg et iend
      // Computation of the Wilkinson's shift
      Real aux = (eigenValues_[iend-1] - eigenValues_[iend])/(2.0 * F_[iend-1]);
      // initialisation of the matrix of rotation
      // y is the current F_[k-1],
      Real y = eigenValues_[ibeg]-eigenValues_[iend] + F_[iend-1]/sign(aux, STK::norm(aux,1.0));
      // z is the element to annulate
      Real z       = F_[ibeg];
      // Fk is the temporary F_[k]
      Real Fk      = z;
      // temporary DeltaD(k)
      Real DeltaDk = 0.0;
      // Index of the columns
      int k,k1;
      // Givens rotation to restore tridiaonal form
      for (k=ibeg, k1=ibeg+1; k<iend; k++, k1++)
      {
        // underflow ? we have a tridiagonal form exit now
        if (z == 0.0) { F_[k]=Fk;  break;}
        // Rotation columns (k,k+1)
        F_[k-1] = (aux = STK::norm(y,z));    // compute final F_[k-1]
        // compute cosinus and sinus
        Real cosinus = y/aux, sinus = -z/aux;
        Real Dk   = eigenValues_[k] + DeltaDk;      // compute current D[k]
        Real aux1 = 2.0 * cosinus * Fk + sinus * (Dk - eigenValues_[k1]);
        // compute DeltaD(k+1)
        eigenValues_[k]     = Dk - (DeltaDk =  sinus * aux1);  // compute eigenValues_[k]
        y         = cosinus*aux1 - Fk;    // temporary F_[k]
        Fk        = cosinus * F_[k1];     // temporary F_[k+1]
        z         = -sinus * F_[k1];      // temporary z
        // update eigenVectors_
        rightGivens(eigenVectors_, k1, k, cosinus, sinus);
      }
      // k=iend if z!=0 and k<iend if z==0
      eigenValues_[k] += DeltaDk ; F_[k-1] = y;
      // restore F[ibeg-1]
      F_[ibeg-1] = 0.;
    } // iter
    if (iter == MAXITER)
    { this->msg_error_ = _T("Warning, max iteration reached in SymEigen::diagonalize()\n");}
    // We have to sort the eigenvalues : we use a basic strategy
    Real z = eigenValues_[iend];        // current value
    for (int i=iend+1; i<eigenValues_.end(); i++)
    { if (eigenValues_[i] > z)           // if the ith eigenvalue is greater
      { eigenValues_.swap(i-1, i);       // swap the cols
        eigenVectors_.swapCols(i-1, i);
      }
      else break;
    } // end sort
  } // iend
  // sort first value
  Real z = eigenValues_[range_.begin()];        // current value
  for (int i=range_.begin()+1; i<eigenValues_.end(); i++)
  { if (eigenValues_[i] > z)           // if the ith eigenvalue is greater
    {
      eigenValues_.swap(i-1, i);       // swap the cols
      eigenVectors_.swapCols(i-1, i);
    }
    else break;
  } // end sort
  // clear auxiliary array
  F_.clear();
}

} // namespace STK

#undef MAXITER

#endif //STK_SYMEIGEN_H
