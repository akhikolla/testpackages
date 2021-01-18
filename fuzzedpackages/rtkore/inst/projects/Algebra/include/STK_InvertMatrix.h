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
 * created on: 9 août 2016
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_InvertMatrix.h
 *  @brief In this file we implement inversion method for general matrix.
 **/

#ifndef STK_INVERTMATRIX_H
#define STK_INVERTMATRIX_H

#include "invert/STK_Invert.h"
#include "invert/STK_InvertUpperTriangular.h"
#include "invert/STK_InvertLowerTriangular.h"
#include "invert/STK_InvertSym.h"
#include "invert/STK_InvertLowerSym.h"
#include "invert/STK_InvertUpperSym.h"

#include <Arrays/include/STK_CArray.h>

namespace STK
{
// forward declaration
template<class Matrix, int Size_> class InvertMatrix;

namespace hidden
{
/** @ingroup hidden
 *  Utility class allowing to choose the dispatcher for inversion matrix
 */
template<class Matrix, int Size_, int structure_>
struct DispatcherChooser;

/** @ingroup hidden
 *  Specialization for general arrays
 */
template<class Matrix, int Size_>
struct DispatcherChooser<Matrix, Size_, Arrays::array2D_>
{
  typedef typename Matrix::Type Type;
  typedef hidden::InvertMatrixDispatcher<Matrix, Size_> Dispatcher;
  typedef CArraySquare<Type, Size_> Result;
};
/** @ingroup hidden
 *  Specialization for square arrays
 */
template<class Matrix, int Size_>
struct DispatcherChooser<Matrix, Size_, Arrays::square_>
{
  typedef typename Matrix::Type Type;
  typedef hidden::InvertMatrixDispatcher<Matrix, Size_> Dispatcher;
  typedef CArraySquare<Type, Size_> Result;
};
/** @ingroup hidden
 *  Specialization for upper triangular arrays
 */
template<class Matrix, int Size_>
struct DispatcherChooser<Matrix, Size_, Arrays::upper_triangular_>
{
  typedef typename Matrix::Type Type;
  typedef hidden::InvertUpperTriangularDispatcher<Matrix, Size_> Dispatcher;
  typedef Array2DUpperTriangular<Type> Result;
};
/** @ingroup hidden
 *  Specialization for upper triangular arrays
 */
template<class Matrix, int Size_>
struct DispatcherChooser<Matrix, Size_, Arrays::lower_triangular_>
{
    typedef typename Matrix::Type Type;
    typedef hidden::InvertLowerTriangularDispatcher<Matrix, Size_> Dispatcher;
    typedef Array2DLowerTriangular<Type> Result;
};
/** @ingroup hidden
 *  Specialization for symmetric arrays
 */
template<class Matrix, int Size_>
struct DispatcherChooser<Matrix, Size_, Arrays::symmetric_>
{
    typedef typename Matrix::Type Type;
    typedef hidden::InvertSymMatrixDispatcher<Matrix, Size_> Dispatcher;
    typedef CArraySquare<Type, Size_> Result;
};
/** @ingroup hidden
 *  Specialization for symmetric arrays
 */
template<class Matrix, int Size_>
struct DispatcherChooser<Matrix, Size_, Arrays::upper_symmetric_>
{
    typedef typename Matrix::Type Type;
    typedef hidden::InvertUpperSymMatrixDispatcher<Matrix, Size_> Dispatcher;
    typedef CArraySquare<Type, Size_> Result;
};
template<class Matrix, int Size_>
struct DispatcherChooser<Matrix, Size_, Arrays::lower_symmetric_>
{
    typedef typename Matrix::Type Type;
    typedef hidden::InvertLowerSymMatrixDispatcher<Matrix, Size_> Dispatcher;
    typedef CArraySquare<Type, Size_> Result;
};

/** @ingroup hidden
 *  Specialization for matrix inversion methods
 */
template<class Matrix, int Size_>
struct AlgebraTraits< InvertMatrix<Matrix, Size_> >
{
  enum
  {
    sizeRows_ = hidden::Traits<Matrix>::sizeRows_,
    sizeCols_ = hidden::Traits<Matrix>::sizeCols_,
    structure_ = hidden::Traits<Matrix>::structure_,
    size_ = (sizeRows_<sizeCols_) ? sizeRows_ : sizeCols_
  };
  typedef typename DispatcherChooser<Matrix, Size_, structure_>::Dispatcher Dispatcher;
  typedef typename DispatcherChooser<Matrix, Size_, structure_>::Result Result;
};


} // namespace hidden

/** @ingroup Algebra
 *  @brief The InvertMatrix class is a functor class allowing to compute the
 *  inverse of a symmetric matrix.
 *
 *  @warning To use with matrices whom first index start at zero with fixed size matrix
 *  from 1 to 4, arbitrary first index is not implemented yet for these matrices.
 */
template<class Matrix, int Size_>
class InvertMatrix
{
  public:
    enum
    {
      size_ = Size_,
      structure_ = hidden::AlgebraTraits< InvertMatrix<Matrix, Size_> >::structure_
    };
    typedef typename hidden::AlgebraTraits<InvertMatrix<Matrix, Size_> >::Dispatcher Dispatcher;
    typedef typename hidden::AlgebraTraits<InvertMatrix<Matrix, Size_> >::Result Result;
    typedef typename Matrix::Type Type;
    /** Constructor */
    inline InvertMatrix( Matrix const& m)
                       : inv_()
                       , det_(Dispatcher::run(m, inv_))
    { inv_.shift(m.beginRows());}
    /** Destructor */
    inline ~InvertMatrix() {}

    /** @return the inverse of the matrix m_ */
    inline Result const& inv() const { return inv_;}
    /** @return the determinant of the matrix m_ */
    inline Type const& det() const { return det_;}
    /** @return @c true if the matrix m_ is invertible, @c false otherwise */
    inline bool isInvertible() const { return(det_!=0);}

    /** compute the inverse of the matrix m_. */
    inline Result const& operator()() { return inv_;}

  protected:
    /** The inverse (or adjugate matrix if det_ is zero) of m_ */
    Result inv_;
    /** determinant of the matrix m_ */
    Type det_;
};

/** @ingroup Algebra
 *  @brief Utility function allowing to compute the inverse of a matrix
 **/
template<class Matrix>
typename hidden::AlgebraTraits< InvertMatrix<Matrix, hidden::Traits<Matrix>::sizeRows_> >::Result
    invert(Matrix const& mat)
{
  enum
  {
    sizeRows_ = hidden::Traits<Matrix>::sizeRows_,
    sizeCols_ = hidden::Traits<Matrix>::sizeCols_,
    size_ = (sizeRows_<sizeCols_) ? sizeRows_ : sizeCols_
  };
  return InvertMatrix<Matrix, sizeRows_>(mat)();
}


}

#endif /* STK_INVERTMATRIX_H_ */
