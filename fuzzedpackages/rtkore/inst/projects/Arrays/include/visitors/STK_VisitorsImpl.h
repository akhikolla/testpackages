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
 * Project:  stkpp::Arrays
 * created on: 27 sept. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_VisitorsImpl.h
 *  @brief In this file we implement the different visitors allowing unrolling
 *  of the visit
 **/

#ifndef STK_VISITORSIMPL_H
#define STK_VISITORSIMPL_H

#define Idx(size) baseIdx + size - 1

namespace STK
{

namespace hidden
{
template<typename Visitor, typename Derived, bool Orient_, int SizeRows_, int SizeCols_>
struct VisitorArrayNoUnrollImpl;
template<typename Visitor, typename Derived, bool Orient_, int SizeRows_, int SizeCols_>
struct VisitorArrayUnrollImpl;
template<typename Visitor, typename Derived, int SizeRows_, int SizeCols_>
struct VisitorArrayImpl;

template<typename Visitor, typename Derived, int SizeRows_, int SizeCols_>
struct VisitorUnrollCol;
template<typename Visitor, typename Derived, int SizeRows_, int SizeCols_>
struct VisitorUnrollRow;

// other kind of structure
template<typename Visitor, typename Derived, int SizeRows_>
struct VisitorVectorImpl;
template<typename Visitor, typename Derived, int SizeRows_>
struct VisitorPointImpl;
template<typename Visitor, typename Derived, int SizeRows_>
struct VisitorDiagonalImpl;
template<typename Visitor, typename Derived, int orient_>
struct VisitorUpperImpl;
template<typename Visitor, typename Derived, int orient_>
struct VisitorLowerImpl;


/** @ingroup hidden
 *  @brief Specialization for general 2D arrays, data stored by column and
 *  dimensions not known at compile time. */
template<typename Visitor, typename Derived>
struct VisitorArrayNoUnrollImpl<Visitor, Derived, Arrays::by_col_, UnknownSize, UnknownSize>
{
  static void run( Derived const& tab, Visitor& visitor)
  {
    for(int j = tab.beginCols(); j < tab.endCols(); ++j)
      for(int i = tab.beginRows(); i < tab.endRows(); ++i)
        visitor(tab.elt(i, j), i, j);
  }
  static void apply( Derived& tab, Visitor& applier)
  {
    for(int j = tab.beginCols(); j < tab.endCols(); ++j)
      for(int i = tab.beginRows(); i < tab.endRows(); ++i)
        applier(tab.elt(i, j));
  }
};

/** @ingroup hidden
 *  @brief Specialization for general 2D arrays, data stored by rows and
 *  dimensions are not known at compile time.*/
template<typename Visitor, typename Derived>
struct VisitorArrayNoUnrollImpl<Visitor, Derived, Arrays::by_row_, UnknownSize, UnknownSize>
{
  static void run( Derived const& tab, Visitor& visitor)
  {
    for(int i = tab.beginRows(); i < tab.endRows(); ++i)
      for(int j = tab.beginCols(); j < tab.endCols(); ++j)
        visitor(tab.elt(i, j), i, j);
  }
  static void apply( Derived& tab, Visitor& applier)
  {
    for(int i = tab.beginRows(); i < tab.endRows(); ++i)
      for(int j = tab.beginCols(); j < tab.endCols(); ++j)
        applier(tab.elt(i, j));
  }
};

/** @ingroup hidden
 *  @brief specialization for the general case when we unroll the rows and the
 *  and the columns with a column oriented arrays */
template<typename Visitor, typename Derived, int SizeRows_ , int SizeCols_>
struct VisitorArrayUnrollImpl<Visitor, Derived, Arrays::by_col_, SizeRows_, SizeCols_>
{
  inline static void run( Derived const& tab, Visitor& visitor)
  {
    // this will unroll the column SizeCols_
    VisitorUnrollCol<Visitor, Derived, SizeRows_, SizeCols_>::run(tab, visitor);
    // do the same for the next column
    VisitorArrayUnrollImpl<Visitor, Derived, Arrays::by_col_, SizeRows_, SizeCols_-1>::run(tab, visitor);
  }
  inline static void apply( Derived& tab, Visitor& applier)
  {
    VisitorUnrollCol<Visitor, Derived, SizeRows_, SizeCols_>::apply(tab, applier);
    VisitorArrayUnrollImpl<Visitor, Derived, Arrays::by_col_, SizeRows_, SizeCols_-1>::apply(tab, applier);
  }
};

/** @ingroup hidden
 *  @brief specialization for the general case when we unroll the rows and the
 *  columns with a row oriented arrays */
template<typename Visitor, typename Derived, int SizeRows_ , int SizeCols_>
struct VisitorArrayUnrollImpl<Visitor, Derived, Arrays::by_row_, SizeRows_, SizeCols_>
{
  inline static void run( Derived const& tab, Visitor& visitor)
  {
    // this will unroll the current row
    VisitorUnrollRow<Visitor, Derived, SizeRows_, SizeCols_>::run(tab, visitor);
    // this will unroll the current row
    VisitorArrayUnrollImpl<Visitor, Derived, Arrays::by_row_, SizeRows_-1, SizeCols_>::run(tab, visitor);
  }
  inline static void apply( Derived& tab, Visitor& visitor)
  {
    VisitorUnrollRow<Visitor, Derived, SizeRows_, SizeCols_>::apply(tab, visitor);
    VisitorArrayUnrollImpl<Visitor, Derived, Arrays::by_row_, SizeRows_-1, SizeCols_>::apply(tab, visitor);
  }
};

/** @ingroup hidden
 *  @brief specialization for the Arrays with 1 column.
 *  In this case, just unroll the column. */
template<typename Visitor, typename Derived, int SizeRows_>
struct VisitorArrayUnrollImpl<Visitor, Derived, Arrays::by_col_, SizeRows_, 1>
{
  inline static void run( Derived const& tab, Visitor& visitor)
  { VisitorUnrollCol<Visitor, Derived, SizeRows_, 1>::run(tab, visitor);}
  inline static void apply( Derived& tab, Visitor& applier)
  { VisitorUnrollCol<Visitor, Derived, SizeRows_, 1>::apply(tab, applier);}
};

/** @ingroup hidden
 *  @brief specialization for the arrays with 1 column (Vector) */
template<typename Visitor, typename Derived, int SizeCols_>
struct VisitorArrayUnrollImpl<Visitor, Derived, Arrays::by_row_, 1, SizeCols_>
{
  inline static void run( Derived const& tab, Visitor& visitor)
  { VisitorUnrollRow<Visitor, Derived, 1, SizeCols_>::run(tab, visitor);}
  inline static void apply( Derived& tab, Visitor& visitor)
  { VisitorUnrollRow<Visitor, Derived, 1, SizeCols_>::apply(tab, visitor);}
};

/** @ingroup hidden
 *  @brief specialization for the general case with 1 row and 1 column arrays.
 *  This specialization allow to disambiguate instantiation.
 **/
template<typename Visitor, typename Derived>
struct VisitorArrayUnrollImpl<Visitor, Derived, Arrays::by_col_, 1, 1>
{
  inline static void run( Derived const& tab, Visitor& visitor)
  { visitor(tab.elt(tab.beginRows(), tab.beginCols()), tab.beginRows(), tab.beginCols());}
  inline static void apply( Derived& tab, Visitor& visitor)
  { visitor(tab.elt(tab.beginRows(), tab.beginCols()));}
};

/** @ingroup hidden
 *  @brief specialization for the general case with 1 row and 1 column arrays.
 *  This specialization allow to disambiguate instantiation.
 **/
template<typename Visitor, typename Derived>
struct VisitorArrayUnrollImpl<Visitor, Derived, Arrays::by_row_, 1, 1>
{
  inline static void run( Derived const& tab, Visitor& visitor)
  { visitor(tab.elt(tab.beginRows(), tab.beginCols()), tab.beginRows(), tab.beginCols());}
  inline static void apply( Derived& tab, Visitor& visitor)
  { visitor(tab.elt(tab.beginRows(), tab.beginCols()));}
};

/** @ingroup hidden
 *  Specialization of VisitorArrayImpl when the number of rows is less than
 *  MaxUnrollSlice and the number of column is unknown or greater
 *  than MaxUnrollSlice
 */
template<typename Visitor, typename Derived, int SizeRows_>
struct VisitorArrayImpl< Visitor, Derived, SizeRows_, UnknownSize>
{
  static void run( Derived const& tab, Visitor& visitor)
  {
    VisitorArrayImpl<Visitor, Derived, SizeRows_ -1, UnknownSize>::run(tab, visitor);
    for(int j = tab.beginCols(); j < tab.endCols(); ++j)
      visitor(tab.elt(Idx(SizeRows_), j), Idx(SizeRows_), j);
  }
  static void apply( Derived& tab, Visitor& visitor)
  {
    VisitorArrayImpl<Visitor, Derived, SizeRows_ -1, UnknownSize>::apply(tab, visitor);
    for(int j = tab.beginCols(); j < tab.endCols(); ++j)
      visitor(tab.elt(Idx(SizeRows_), j));
  }
};

/** @ingroup hidden
 *  @brief Specialization of VisitorArrayImpl when the number of rows is 1
 *  and the number of column is unknown or greater than MaxUnrollSlice.
 *  This class will stop the unrolling of on the rows.
 **/
template<typename Visitor, typename Derived>
struct VisitorArrayImpl< Visitor, Derived, 1, UnknownSize>
{
  static void run( Derived const& tab, Visitor& visitor)
  {
    for(int j = tab.beginCols(); j < tab.endCols(); ++j)
      visitor(tab.elt(tab.beginRows(), j), tab.beginRows(), j);
  }
  static void apply( Derived& tab, Visitor& visitor)
  {
    for(int j = tab.beginCols(); j < tab.endCols(); ++j)
      visitor(tab.elt(tab.beginRows(), j));
  }
};

/** @ingroup hidden
 *  Specialization of VisitorArrayImpl when the number of columns is less than
 *  MaxUnrollSlice and the number of rows is unknown
 */
template<typename Visitor, typename Derived, int SizeCols_>
struct VisitorArrayImpl< Visitor, Derived, UnknownSize, SizeCols_>
{
  static void run( Derived const& tab, Visitor& visitor)
  {
    VisitorArrayImpl<Visitor, Derived, UnknownSize, SizeCols_-1>::run(tab, visitor);
    for(int i = tab.beginRows(); i < tab.endRows(); ++i)
      visitor(tab.elt(i, Idx(SizeCols_)), i, Idx(SizeCols_));
  }
  static void apply( Derived& tab, Visitor& visitor)
  {
    VisitorArrayImpl<Visitor, Derived, UnknownSize, SizeCols_ -1>::apply(tab, visitor);
    for(int i = tab.beginRows(); i < tab.endRows(); ++i)
      visitor(tab.elt(i, Idx(SizeCols_)));
  }
};

/** @ingroup hidden
 *  Specialization of VisitorArrayImpl when the number of columns is one and
 *  the number of rows is unknown.
 */
template<typename Visitor, typename Derived>
struct VisitorArrayImpl< Visitor, Derived, UnknownSize, 1>
{
  static void run( Derived const& tab, Visitor& visitor)
  {
    for(int i = tab.beginRows(); i < tab.endRows(); ++i)
      visitor(tab.elt(i, tab.beginCols()), i, tab.beginCols());
  }
  static void apply( Derived& tab, Visitor& applier)
  {
    for(int i = tab.beginRows(); i < tab.endRows(); ++i)
      applier(tab.elt(i, tab.beginCols()));
  }
};

/** @ingroup hidden
 *  @brief unroll the column baseIdx + SizeCols_ -1 */
template<typename Visitor, typename Derived, int SizeRows_, int SizeCols_>
struct VisitorUnrollCol
{
  inline static void run( Derived const& tab, Visitor& visitor)
  {
    VisitorUnrollCol<Visitor, Derived, SizeRows_-1, SizeCols_>::run(tab, visitor);
    visitor(tab.elt(Idx(SizeRows_), Idx(SizeCols_)), Idx(SizeRows_), Idx(SizeCols_));
  }
  inline static void apply( Derived& tab, Visitor& visitor)
  {
    VisitorUnrollCol<Visitor, Derived, SizeRows_-1, SizeCols_>::apply(tab, visitor);
    visitor(tab.elt(Idx(SizeRows_), Idx(SizeCols_)));
  }
};

/** @ingroup hidden
 *  @brief specialization for the arrays with 1 row (Point) */
template<typename Visitor, typename Derived, int SizeCols_>
struct VisitorUnrollCol<Visitor, Derived, 1, SizeCols_>
{
  inline static void run( Derived const& tab, Visitor& visitor)
  { visitor(tab.elt(baseIdx, Idx(SizeCols_)), baseIdx, Idx(SizeCols_));}
  inline static void apply( Derived& tab, Visitor& visitor)
  { visitor(tab.elt(baseIdx, Idx(SizeCols_)));}
};

/** @ingroup hidden
 *  @brief  unroll the row baseIdx + SizeRows_ -1  */
template<typename Visitor, typename Derived, int SizeRows_, int SizeCols_>
struct VisitorUnrollRow
{
  inline static void run( Derived const& tab, Visitor& visitor)
  {
    VisitorUnrollRow<Visitor, Derived, SizeRows_, SizeCols_-1>::run(tab, visitor);
    visitor(tab.elt(Idx(SizeRows_), Idx(SizeCols_)), Idx(SizeRows_), Idx(SizeCols_));
  }
  inline static void apply( Derived& tab, Visitor& visitor)
  {
    VisitorUnrollRow<Visitor, Derived, SizeRows_, SizeCols_-1>::apply(tab, visitor);
    visitor(tab.elt(Idx(SizeRows_), Idx(SizeCols_)));
  }
};

/** @ingroup hidden
 *  @brief specialization for the arrays with 1 column */
template<typename Visitor, typename Derived, int SizeRows_>
struct VisitorUnrollRow<Visitor, Derived, SizeRows_, 1>
{
  inline static void run( Derived const& tab, Visitor& visitor)
  { visitor(tab.elt(Idx(SizeRows_), tab.beginCols()), Idx(SizeRows_), tab.beginCols());}
  inline static void apply( Derived& tab, Visitor& visitor)
  { visitor(tab.elt(Idx(SizeRows_), tab.beginCols()));}
};

/** @ingroup hidden
 *  @brief Specialization when the size is unknown */
template<typename Visitor, typename Derived>
struct VisitorVectorImpl<Visitor, Derived, UnknownSize>
{
  static void run( Derived const& tab, Visitor& visitor)
  { for(int i = tab.begin(); i < tab.end(); ++i)
    visitor(tab.elt(i), i, tab.colIdx());}
  static void apply( Derived& tab, Visitor& applier)
  { for(int i = tab.begin(); i < tab.end(); ++i) applier(tab.elt(i));}
};

/** @ingroup hidden
 *  @brief A visitor Vector allow to unroll the visit of a vector if possible
 **/
template<typename Visitor, typename Derived, int Size_>
struct VisitorVectorImpl
{
  inline static void run( Derived const& tab, Visitor& visitor)
  {
    VisitorVectorImpl<Visitor, Derived, Size_-1>::run(tab, visitor);
    visitor(tab.elt(tab.begin()+ Size_-1), tab.begin()+ Size_-1, tab.colIdx());
  }
  inline static void apply( Derived& tab, Visitor& applier)
  {
    VisitorVectorImpl<Visitor, Derived, Size_-1>::apply(tab, applier);
    applier(tab.elt(tab.begin()+ Size_-1));
  }
};

/** @ingroup hidden
 *  @brief specialization for the vectors with 1 row */
template<typename Visitor, typename Derived>
struct VisitorVectorImpl<Visitor, Derived, 1>
{
  inline static void run( Derived const& tab, Visitor& visitor)
  { visitor(tab.elt(tab.begin()), tab.begin(), tab.colIdx());}
  inline static void apply( Derived& tab, Visitor& applier)
  { applier(tab.elt(tab.begin()));}
};

/** @ingroup hidden
 *  @brief Specialization when the size is unknown */
template<typename Visitor, typename Derived>
struct VisitorPointImpl<Visitor, Derived, UnknownSize>
{
  static void run( Derived const& tab, Visitor& visitor)
  { for(int j = tab.begin(); j < tab.end(); ++j) visitor(tab.elt(j), tab.rowIdx(),j);}
  static void apply( Derived& tab, Visitor& applier)
  { for(int j = tab.begin(); j < tab.end(); ++j) applier(tab.elt(j));}
};


/** @ingroup hidden
 *  @brief A VisitorPointImpl allow to unroll the visit of a row-vector if possible
 *  */
template<typename Visitor, typename Derived, int Size_>
struct VisitorPointImpl
{
  inline static void run( Derived const& tab, Visitor& visitor)
  {
    VisitorPointImpl<Visitor, Derived, Size_-1>::run(tab, visitor);
    visitor(tab.elt(tab.begin()+ Size_-1), tab.rowIdx(), tab.begin()+ Size_-1);
  }
  inline static void apply( Derived& tab, Visitor& visitor)
  {
    VisitorPointImpl<Visitor, Derived, Size_-1>::apply(tab, visitor);
    visitor(tab.elt(tab.begin()+ Size_-1));
  }
};

/** @ingroup hidden
 *  @brief specialization for the RowVectors with one columns */
template<typename Visitor, typename Derived>
struct VisitorPointImpl<Visitor, Derived, 1>
{
  inline static void run( Derived const& tab, Visitor& visitor)
  { visitor(tab.elt(tab.begin()), tab.rowIdx(), tab.begin());}
  inline static void apply( Derived& tab, Visitor& visitor)
  { visitor(tab.elt(tab.begin()));}
};

/** @ingroup hidden
 *  @brief Specialization when the size is unknown */
template<typename Visitor, typename Derived>
struct VisitorDiagonalImpl<Visitor, Derived, UnknownSize>
{
  static void run( Derived const& tab, Visitor& visitor)
  { for(int j = tab.begin(); j < tab.end(); ++j) visitor(tab.elt(j), j,j);}
  static void apply( Derived& tab, Visitor& visitor)
  { for(int j = tab.begin(); j < tab.end(); ++j) visitor(tab.elt(j));}
};

/** @ingroup hidden
 *  @brief A VisitorDiagonalImpl allow to unroll the visit of a Diagonal tab if possible
 *  */
template<typename Visitor, typename Derived, int Size_>
struct VisitorDiagonalImpl
{
  inline static void run( Derived const& tab, Visitor& visitor)
  {
    VisitorPointImpl<Visitor, Derived, Size_-1>::run(tab, visitor);
    visitor(tab.elt(tab.begin()+ Size_-1), tab.begin()+ Size_-1, tab.begin()+ Size_-1);
  }
  inline static void apply( Derived& tab, Visitor& visitor)
  {
    VisitorPointImpl<Visitor, Derived, Size_-1>::apply(tab, visitor);
    visitor(tab.elt(tab.begin()+ Size_-1));
  }
};

/** @ingroup hidden
 *  @brief specialization for diagonal tab of size one */
template<typename Visitor, typename Derived>
struct VisitorDiagonalImpl<Visitor, Derived, 1>
{
  inline static void run( Derived const& tab, Visitor& visitor)
  { visitor(tab.elt(tab.begin()), tab.begin(), tab.begin());}
  inline static void apply( Derived& tab, Visitor& visitor)
  { visitor(tab.elt(tab.begin()));}
};

/** @ingroup hidden
 *  @brief specialization for the general case with column oriented arrays */
template<typename Visitor, typename Derived>
struct VisitorUpperImpl<Visitor, Derived, Arrays::by_col_>
{
  static void run( Derived const& tab, Visitor& visitor)
  {
    for(int j = tab.lastIdxCols(); j >= tab.beginCols(); --j)
      for(int i = std::min(j, tab.lastIdxRows()); i >= tab.beginRows(); --i)
        visitor(tab.elt(i, j), i, j);
  }
  static void apply( Derived& tab, Visitor& visitor)
  {
    for(int j = tab.lastIdxCols(); j >= tab.beginCols(); --j)
      for(int i = std::min(j, tab.lastIdxRows()); i >= tab.beginRows(); --i)
        visitor(tab.elt(i, j));
  }
};

/** @ingroup hidden
 *  @brief specialization for the general case with row oriented arrays */
template<typename Visitor, typename Derived>
struct VisitorUpperImpl<Visitor, Derived, Arrays::by_row_>
{
  static void run( Derived const& tab, Visitor& visitor)
  {
    for(int i = tab.beginRows(); i < tab.endRows(); ++i)
      for(int j = std::max(i, tab.beginCols()); j < tab.endCols(); ++j)
        visitor(tab.elt(i, j), i, j);
  }
  static void apply( Derived& tab, Visitor& visitor)
  {
    for(int i = tab.beginRows(); i < tab.endRows(); ++i)
      for(int j = std::max(i, tab.beginCols()); j < tab.endCols(); ++j)
        visitor(tab.elt(i, j));
  }
};


/** @ingroup hidden
 *  @brief specialization for the general case with column oriented arrays */
template<typename Visitor, typename Derived>
struct VisitorLowerImpl<Visitor, Derived, Arrays::by_col_>
{
  static void run( Derived const& tab, Visitor& visitor)
  {
    for(int j = tab.beginCols(); j < tab.endCols(); ++j)
      for(int i = std::max(j, tab.beginRows()); i < tab.endRows(); ++i)
        visitor(tab.elt(i, j), i, j);
  }
  static void apply( Derived& tab, Visitor& visitor)
  {
    for(int j = tab.beginCols(); j < tab.endCols(); ++j)
      for(int i = std::max(j, tab.beginRows()); i < tab.endRows(); ++i)
        visitor(tab.elt(i, j));
  }
};

/** @ingroup hidden
 *  @brief specialization for the general case with row oriented arrays */
template<typename Visitor, typename Derived>
struct VisitorLowerImpl<Visitor, Derived, Arrays::by_row_>
{
  static void run( Derived const& tab, Visitor& visitor)
  {
    for(int i = tab.lastIdxRows(); i >= tab.beginRows(); --i)
      for(int j = std::min(i, tab.lastIdxCols()); j >= tab.beginCols(); --j)
        visitor(tab.elt(i, j), i, j);
  }
  static void apply( Derived& tab, Visitor& visitor)
  {
    for(int i = tab.lastIdxRows(); i >= tab.beginRows(); --i)
      for(int j = std::min(i, tab.lastIdxCols()); j >= tab.beginCols(); --j)
        visitor(tab.elt(i, j));
  }
};

/** @ingroup hidden
 *  @brief specialization for the general case with 1 row and 1 column arrays.
 *  This specialization allow to disambiguate instantiation.
 **/
template<typename Visitor, typename Derived>
struct VisitorNumberImpl
{
  inline static void run( Derived const& tab, Visitor& visitor)
  { visitor(tab.elt(), tab.beginRows(), tab.beginCols());}
  inline static void apply( Derived& tab, Visitor& visitor)
  { visitor(tab.elt());}
};


} // namespace hidden

} // namespace STK

#undef Idx

#endif /* STK_VISITORSIMPL_H */
