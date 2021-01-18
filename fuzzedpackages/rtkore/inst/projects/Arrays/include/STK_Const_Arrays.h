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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Arrays
 * created on: 25 juil. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Const_Arrays.h
 *  @brief In this file we define the constant Arrays.
 **/


#ifndef STK_CONST_ARRAYS_H
#define STK_CONST_ARRAYS_H

#include "STK_ExprBase.h"

namespace STK
{

// forward declaration
namespace Const
{
template< typename Type_, int Size_ = UnknownSize> class Identity;
template< typename Type_, int Size_ = UnknownSize> class Square;
template< typename Type_, int Size_ = UnknownSize> class Vector;
template< typename Type_, int Size_ = UnknownSize> class Point;
template< typename Type_, int SizeRows_ = UnknownSize, int SizeCols_ = UnknownSize> class Array;
template< typename Type_, int SizeRows_ = UnknownSize, int SizeCols_ = UnknownSize> class UpperTriangular;
template< typename Type_, int SizeRows_ = UnknownSize, int SizeCols_ = UnknownSize> class LowerTriangular;

typedef Vector<Real, UnknownSize>   VectorX;
typedef Vector<Real, 3>             Vector3;
typedef Vector<Real, 2>             Vector2;
typedef Vector<double, UnknownSize> VectorXd;
typedef Vector<double, 3>           Vector3d;
typedef Vector<double, 2>           Vector2d;
typedef Vector<int, UnknownSize>    VectorXi;
typedef Vector<int, 3>              Vector3i;
typedef Vector<int, 2>              Vector2i;

typedef Vector<Real, UnknownSize>   ConstVectorX;
typedef Vector<Real, 3>             ConstVector3;
typedef Vector<Real, 2>             ConstVector2;
typedef Vector<double, UnknownSize> ConstVectorXd;
typedef Vector<double, 3>           ConstVector3d;
typedef Vector<double, 2>           ConstVector2d;
typedef Vector<int, UnknownSize>    ConstVectorXi;
typedef Vector<int, 3>              ConstVector3i;
typedef Vector<int, 2>              ConstVector2i;

typedef Point<Real, UnknownSize>   PointX;
typedef Point<Real, 3>             Point3;
typedef Point<Real, 2>             Point2;
typedef Point<double, UnknownSize> PointXd;
typedef Point<double, 3>           Point3d;
typedef Point<double, 2>           Point2d;
typedef Point<int, UnknownSize>    PointXi;
typedef Point<int, 3>              Point3i;
typedef Point<int, 2>              Point2i;

typedef Point<Real, UnknownSize>   ConstPointX;
typedef Point<Real, 3>             ConstPoint3;
typedef Point<Real, 2>             ConstPoint2;
typedef Point<double, UnknownSize> ConstPointXd;
typedef Point<double, 3>           ConstPoint3d;
typedef Point<double, 2>           ConstPoint2d;
typedef Point<int, UnknownSize>    ConstPointXi;
typedef Point<int, 3>              ConstPoint3i;
typedef Point<int, 2>              ConstPoint2i;
}

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the identity constant matrix
 */
template< typename Type_, int Size_>
struct Traits< Const::Identity <Type_, Size_> >
{
  enum
  {
    structure_ = Arrays::diagonal_,
    sizeRows_  = Size_,
    sizeCols_  = Size_,
    size_      = Size_,
    orient_    = Arrays::by_col_,
    storage_   = Arrays::sparse_
  };

  typedef Type_                Type;
  typedef typename RemoveConst<Type>::Type ConstReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for the identity constant matrix
 */
template< typename Type_, int Size_>
struct Traits< Const::Vector <Type_, Size_> >
{
  enum
  {
    structure_ = Arrays::vector_,
    sizeRows_  = Size_,
    sizeCols_  = 1,
    size_      = Size_,
    orient_    = Arrays::by_col_,
    storage_   = Arrays::sparse_
  };
  typedef Type_                Type;
  typedef typename RemoveConst<Type>::Type ConstReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for the constant vector_
 */
template< typename Type_, int Size_>
struct Traits< Const::Point <Type_, Size_> >
{
  enum
  {
    structure_ = Arrays::point_,
    sizeRows_  = 1,
    sizeCols_  = Size_,
    size_      = Size_,
    orient_    = Arrays::by_row_,
    storage_   = Arrays::sparse_
  };
  typedef Type_                Type;
  typedef typename RemoveConst<Type>::Type ConstReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for the square constant matrix
 */
template< typename Type_, int Size_>
struct Traits< Const::Square <Type_, Size_> >
{
  enum
  {
    structure_ = Arrays::square_,
    sizeRows_  = Size_,
    sizeCols_  = Size_,
    size_      = Size_,
    orient_    = Arrays::by_col_,
    storage_   = Arrays::sparse_
  };
  typedef Type_                Type;
  typedef typename RemoveConst<Type>::Type ConstReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for the general constant matrix
 */
template< typename Type_, int SizeRows_, int SizeCols_>
struct Traits< Const::Array<Type_, SizeRows_, SizeCols_> >
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_  = SizeRows_,
    sizeCols_  = SizeCols_,
    orient_    = Arrays::by_col_,
    storage_   = Arrays::sparse_
  };
  typedef Type_                Type;
  typedef typename RemoveConst<Type>::Type ConstReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for the upper triangular constant matrix
 */
template< typename Type_, int SizeRows_, int SizeCols_>
struct Traits< Const::UpperTriangular<Type_, SizeRows_, SizeCols_> >
{
  enum
  {
    structure_ = Arrays::upper_triangular_,
    sizeRows_  = SizeRows_,
    sizeCols_  = SizeCols_,
    orient_    = Arrays::by_col_,
    storage_   = Arrays::sparse_
  };
  typedef Type_                Type;
  typedef typename RemoveConst<Type>::Type ConstReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for the lower triangular constant matrix
 */
template< typename Type_, int SizeRows_, int SizeCols_>
struct Traits< Const::LowerTriangular<Type_, SizeRows_, SizeCols_> >
{
  enum
  {
    structure_ = Arrays::lower_triangular_,
    sizeRows_  = SizeRows_,
    sizeCols_  = SizeCols_,
    orient_    = Arrays::by_col_,
    storage_   = Arrays::sparse_
  };
  typedef Type_                Type;
  typedef typename RemoveConst<Type>::Type ConstReturnType;
};

} // namespace hidden

namespace Const
{
template< class Derived>
class IConstArray: protected IContainer2D<hidden::Traits<Derived>::sizeRows_, hidden::Traits<Derived>::sizeCols_>, public ExprBase<Derived>
{
  public:
    enum
    {
      structure_ = hidden::Traits< Derived >::structure_,
      orient_    = hidden::Traits< Derived >::orient_,
      sizeRows_  = hidden::Traits< Derived >::sizeRows_,
      sizeCols_  = hidden::Traits< Derived >::sizeCols_,
      storage_   = hidden::Traits< Derived >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;
    /** Type for the IContainer2D base Class. */
    typedef IContainer2D<sizeRows_, sizeCols_ > Base2D;
    typedef ExprBase<Derived> Base;
    typedef typename hidden::Traits< Derived >::ConstReturnType ConstReturnType;

  protected:
    /** default constructor */
    IConstArray(): Base2D(), Base() {}
    /** constructor with specified dimension */
    IConstArray(Range const& rows, Range const& cols): Base2D(rows, cols), Base() {}
    /** destructor */
    ~IConstArray() {}

  public:
    /** @return the Vertical range */
    inline RowRange const& rowsImpl() const { return Base2D::rows();}
    /** @return the Vertical range */
    inline RowRange const& rows() const { return Base2D::rows();}
    /** @return the index of the first row */
    inline int beginRows() const { return Base2D::beginRows();}
    /** @return the ending index of the rows */
    inline int endRows() const { return Base2D::endRows();}
    /** @return the number of rows */
    inline int sizeRows() const { return Base2D::sizeRows();}

    /**@return the Horizontal range */
    inline ColRange const&colsImpl() const { return Base2D::cols();}
    /**@return the Horizontal range */
    inline ColRange const&cols() const { return Base2D::cols();}
    /** @return the index of the first column */
    inline int beginCols() const { return Base2D::beginCols();}
    /**  @return the ending index of columns */
    inline int endCols() const { return Base2D::endCols();}
    /** @return the number of columns */
    inline int sizeCols() const { return Base2D::sizeCols();}

    /**  @return the index of the last column */
    inline int lastIdxCols() const { return Base2D::lastIdxCols();}
    /** @return the index of the last row */
    inline int lastIdxRows() const { return Base2D::lastIdxRows();}

    /**  @return @c true if the container is empty, @c false otherwise */
    bool empty() const { return Base2D::empty();}
};

/**@ingroup Arrays
 * Define the constant identity matrix
 * \f[
 * I_n =
 * \left(
 * \begin{array}{ccccc}
 *  1\\
 *   & 1 & & \text{\huge0}\\
 *   &   & \ddots \\
 *   & \text{\huge0} & & \ddots\\
 *   &               & &  & 1
 * \end{array}
 * \right).
 * \f]
 * The size can be either a fixed template argument or a dynamic size.
 * Exemple:
 * @code
 *  STK::Square<Real,3> S3; // S3 is a 3x3 square matrix of 1.
 *  STK::Square<Real> S(10); // S is a 10x10 identity matrix of 1.
 * @endcode
 * @tparam Size_ the size of the identity ArrayXX. Default is UnknownSize.
 **/
template< typename Type_, int Size_ >
class Identity: public IConstArray<Identity<Type_, Size_> >
{
  public:
    typedef IConstArray<Identity<Type_, Size_> > Base;
    typedef typename hidden::Traits< Identity<Type_, Size_> >::Type Type;
    typedef typename hidden::Traits< Identity<Type_, Size_> >::ConstReturnType ConstReturnType;
    enum
    {
      structure_ = hidden::Traits< Identity <Type_, Size_> >::structure_,
      orient_    = hidden::Traits< Identity <Type_, Size_> >::orient_,
      sizeRows_  = hidden::Traits< Identity <Type_, Size_> >::sizeRows_,
      sizeCols_  = hidden::Traits< Identity <Type_, Size_> >::sizeCols_,
      storage_   = hidden::Traits< Identity <Type_, Size_> >::storage_
    };
    /** default constructor */
    Identity(): Base() {}
    /** constructor with specified dimension */
    Identity(Range I): Base(I,I) {}
    /** @return the element (i,j) of the identity matrix.
     *  @param i,j row and column indexes
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (i==j ? Type(1) : Type(0));}
    /** @return the element ith element of the identity matrix
     *  @param i index of the element to get
     **/
    inline ConstReturnType elt1Impl(int i) const { return Type(1);}
};

/**@ingroup Arrays
 * Define the constant square matrix
 * \f[
 * S =
 * \left(
 * \begin{array}{ccccc}
 *  1\\
 *   & 1 & & \text{\huge 1}\\
 *   &   & \ddots \\
 *   & \text{\huge 1} & & \ddots\\
 *   &               & &  & 1
 * \end{array}
 * \right).
 * \f]
 * The size can be either a fixed template argument or a dynamic size.
 * Exemple:
 * @code
 *  STK::Square<Real,3> S3; // S3 is a 3x3 square matrix of Real
 *  STK::Square<Real> S(10); // S is a 10x10 square matrix of Real
 * @endcode
 * @tparam Size_ the size of the square ArrayXX. Default is UnknownSize.
 **/
template< typename Type_, int Size_ >
class Square: public IConstArray<Square<Type_, Size_> >
{
  public:
    typedef IConstArray<Square<Type_, Size_> > Base;
    typedef typename hidden::Traits< Square<Type_, Size_> >::Type Type;
    typedef typename hidden::Traits< Square<Type_, Size_> >::ConstReturnType ConstReturnType;
    enum
    {
      structure_ = hidden::Traits< Square <Type_, Size_> >::structure_,
      orient_    = hidden::Traits< Square <Type_, Size_> >::orient_,
      sizeRows_  = hidden::Traits< Square <Type_, Size_> >::sizeRows_,
      sizeCols_  = hidden::Traits< Square <Type_, Size_> >::sizeCols_,
      storage_   = hidden::Traits< Square <Type_, Size_> >::storage_
    };
    /** default constructor */
    Square(): Base() {}
    /** constructor with specified dimension */
    Square(Range I): Base(I, I) {}
    /** @return the element (i,j) of the constant square matrix.
     *  @param i,j row and column indexes
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (Type(1));}
};

/**@ingroup Arrays
 * Define the constant general matrix
 * \f[
 * S =
 * \left(
 * \begin{array}{ccccc}
 *  1      &  & \hdots & & 1\\
 *         &  &        & &  \\
 *  \vdots &  &        & & \vdots\\
 *         &  &        & & \\
 *  1      &  & \hdots & & 1
 * \end{array}
 * \right).
 * \f]
 * The sizes can be either two fixed template arguments or dynamic sizes.
 * Exemple:
 * @code
 *  STK::Array<Real,3, 4> G3; // G3 is a 3x4 matrix of Real
 *  STK::Array<Real> S(10, 20); // S is a 10x20 matrix of Real
 * @endcode
 * @tparam SizeRows_ the number of row of the matrix. Default is UnknownSize.
 * @tparam SizeCols_ the number of column of the matrix. Default is UnknownSize.
 **/
template< typename Type_, int SizeRows_, int SizeCols_ >
class Array: public IConstArray<Array<Type_, SizeRows_, SizeCols_> >
{
  public:
    typedef IConstArray<Array<Type_, SizeRows_, SizeCols_> > Base;
    typedef typename hidden::Traits< Array<Type_, SizeRows_, SizeCols_> >::Type Type;
    typedef typename hidden::Traits< Array<Type_, SizeRows_, SizeCols_> >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< Array <Type_, SizeRows_, SizeCols_> >::structure_,
      orient_    = hidden::Traits< Array <Type_, SizeRows_, SizeCols_> >::orient_,
      sizeRows_  = hidden::Traits< Array <Type_, SizeRows_, SizeCols_> >::sizeRows_,
      sizeCols_  = hidden::Traits< Array <Type_, SizeRows_, SizeCols_> >::sizeCols_,
      storage_   = hidden::Traits< Array <Type_, SizeRows_, SizeCols_> >::storage_
    };
    /** default constructor */
    Array(): Base() {}
    /** constructor with specified dimension */
    Array(Range rangeRows, Range rangeCols): Base(rangeRows, rangeCols) {}
    /** @return the element (i,j) of the constant square matrix.
     *  @param i,j row and column indexes
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (Type(1));}
};

/**@ingroup Arrays
 * Define the constant upper triangular matrix
 * \f[
 * S =
 * \left(
 * \begin{array}{ccccc}
 *  1\\
 *   & 1 & & \text{\huge 1}\\
 *   &   & \ddots \\
 *   & \text{\huge 0} & & \ddots\\
 *   &               & &  & 1
 * \end{array}
 * \right).
 * \f]
 * The sizes can be either two fixed template arguments or dynamic sizes.
 * Exemple:
 * @code
 *  STK::UpperTriangular<Real,3, 4> U3; // U3 is a 3x4 upper triangular matrix of Real
 *  STK::UpperTriangular<Real> S(10, 20); // S is a 10x20 upper triangula matrix of Real
 * @endcode
 * @tparam SizeRows_ the number of row of the matrix. Default is UnknownSize.
 * @tparam SizeCols_ the number of column of the matrix. Default is UnknownSize.
 **/
template< typename Type_, int SizeRows_, int SizeCols_ >
class UpperTriangular: public IConstArray<UpperTriangular<Type_, SizeRows_, SizeCols_> >
{
  public:
    typedef IConstArray<UpperTriangular<Type_, SizeRows_, SizeCols_> > Base;
    typedef typename hidden::Traits< UpperTriangular<Type_, SizeRows_, SizeCols_> >::Type Type;
    typedef typename hidden::Traits< UpperTriangular<Type_, SizeRows_, SizeCols_> >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< UpperTriangular <Type_, SizeRows_, SizeCols_> >::structure_,
      orient_    = hidden::Traits< UpperTriangular <Type_, SizeRows_, SizeCols_> >::orient_,
      sizeRows_  = hidden::Traits< UpperTriangular <Type_, SizeRows_, SizeCols_> >::sizeRows_,
      sizeCols_  = hidden::Traits< UpperTriangular <Type_, SizeRows_, SizeCols_> >::sizeCols_,
      storage_   = hidden::Traits< UpperTriangular <Type_, SizeRows_, SizeCols_> >::storage_
    };
    /** default constructor */
    UpperTriangular(): Base() {}
    /** constructor with specified dimension */
    UpperTriangular(Range rangeRows, Range rangeCols): Base(rangeRows, rangeCols) {}
    /** @return the element (i,j) of the constant upper triangular matrix.
     *  @param i,j row and column indexes
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return i<j ? (Type(1)) : Type(0);}
};

/**@ingroup Arrays
 * Define the constant lower triangular matrix
 * \f[
 * S =
 * \left(
 * \begin{array}{ccccc}
 *  1\\
 *   & 1 & & \text{\huge 0}\\
 *   &   & \ddots \\
 *   & \text{\huge 1} & & \ddots\\
 *   &               & &  & 1
 * \end{array}
 * \right).
 * \f]
 * The sizes can be either two fixed template arguments or dynamic sizes.
 * Exemple:
 * @code
 *  STK::LowerTriangular<Real,3, 4> G3; // G3 is a 3x4 lower triangular matrix of Real
 *  STK::LowerTriangular<Real> S(10, 20); // S is a 10x20 lower triangular matrix of Real
 * @endcode
 * @tparam SizeRows_ the number of row of the matrix. Default is UnknownSize.
 * @tparam SizeCols_ the number of column of the matrix. Default is UnknownSize.
 **/
template< typename Type_, int SizeRows_, int SizeCols_ >
class LowerTriangular: public IConstArray<LowerTriangular<Type_, SizeRows_, SizeCols_> >
{
  public:
    typedef IConstArray<LowerTriangular<Type_, SizeRows_, SizeCols_> > Base;
    typedef typename hidden::Traits< LowerTriangular<Type_, SizeRows_, SizeCols_> >::Type Type;
    typedef typename hidden::Traits< LowerTriangular<Type_, SizeRows_, SizeCols_> >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< LowerTriangular <Type_, SizeRows_, SizeCols_> >::structure_,
      orient_    = hidden::Traits< LowerTriangular <Type_, SizeRows_, SizeCols_> >::orient_,
      sizeRows_  = hidden::Traits< LowerTriangular <Type_, SizeRows_, SizeCols_> >::sizeRows_,
      sizeCols_  = hidden::Traits< LowerTriangular <Type_, SizeRows_, SizeCols_> >::sizeCols_,
      storage_   = hidden::Traits< LowerTriangular <Type_, SizeRows_, SizeCols_> >::storage_
    };
    /** default constructor */
    LowerTriangular(): Base() {}
    /** constructor with specified dimension */
    LowerTriangular(Range rangeRows, Range rangeCols): Base(rangeRows, rangeCols) {}
    /** @return the element (i,j) of the constant lower triangular matrix.
     *  @param i,j row and column indexes
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (Type(1));}
};

/**@ingroup Arrays
 * Define the constant point
 * \f[
 * S =
 * \left(
 * \begin{array}{cccc}
 *  1 & 1 & \ldots  & 1
 * \end{array}
 * \right).
 * \f]
 * The size can be either a fixed template argument or a dynamic size.
 * Exemple:
 * @code
 *  STK::Point<Real,3> p3; // p3 is a row-vector of size 3
 *  STK::Point<Real> p(10); // p is a row-vector of Real of size 10
 * @endcode
 * @tparam Size_ the size of the row-vector. Default is UnknownSize.
 **/
template< typename Type_, int Size_ >
class Point: public IConstArray<Point<Type_, Size_> >
{
  public:
    typedef IConstArray<Point<Type_, Size_> > Base;
    typedef typename hidden::Traits< Point<Type_, Size_> >::Type Type;
    typedef typename hidden::Traits< Point<Type_, Size_> >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< Point <Type_, Size_> >::structure_,
      orient_    = hidden::Traits< Point <Type_, Size_> >::orient_,
      sizeRows_  = hidden::Traits< Point <Type_, Size_> >::sizeRows_,
      sizeCols_  = hidden::Traits< Point <Type_, Size_> >::sizeCols_,
      storage_   = hidden::Traits< Point <Type_, Size_> >::storage_
    };
    /** default constructor */
    Point(): Base() {}
    /** constructor with specified dimension */
    Point(Range I): Base(1, I) {}
    /** @return the j-th element  of the constant row-vector.
     *  @param j index of the element
     **/
    inline ConstReturnType elt1Impl(int j) const { return (Type(1));}
    /** @return the element (i,j) of the constant row-vector.
     *  @param i,j row and column indexes
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (Type(1));}
};

/**@ingroup Arrays
 * Define the constant point
 * \f[
 * S =
 * \left(
 * \begin{array}{cccc}
 *  1 & 1 & \ldots  & 1
 * \end{array}
 * \right).
 * \f]
 * The size can be either a fixed template argument or a dynamic size.
 * Exemple:
 * @code
 *  STK::Point<Real,3> p3; // p3 is a row-vector of size 3
 *  STK::Point<Real> p(10); // p is a row-vector of Real of size 10
 * @endcode
 * @tparam Size_ the size of the row-vector. Default is UnknownSize.
 **/
template< typename Type_, int Size_ >
class Vector: public IConstArray<Vector<Type_, Size_> >
{
  public:
    typedef IConstArray<Vector<Type_, Size_> > Base;
    typedef typename hidden::Traits< Vector<Type_, Size_> >::Type Type;
    typedef typename hidden::Traits< Vector<Type_, Size_> >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< Vector <Type_, Size_> >::structure_,
      orient_    = hidden::Traits< Vector <Type_, Size_> >::orient_,
      sizeRows_  = hidden::Traits< Vector <Type_, Size_> >::sizeRows_,
      sizeCols_  = hidden::Traits< Vector <Type_, Size_> >::sizeCols_,
      storage_   = hidden::Traits< Vector <Type_, Size_> >::storage_
    };
    /** default constructor */
    Vector(): Base() {}
    /** constructor with specified dimension */
    Vector(Range I): Base(I, 1) {}
    /** @return the i-th element of the constant vector.
     *  @param i index of the element
     **/
    inline ConstReturnType elt1Impl(int i) const { return (Type(1));}
    /** @return the element (i,j) of the constant vector.
     *  @param i,j row and column indexes
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (Type(1));}
};


} // namespace const

} // namespace STK

#endif /* STK_CONST_ARRAYS_H */
