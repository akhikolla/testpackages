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
 * Purpose:  Functors applied to the Expression and arrays
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_SlicingVisitors.h
 *  @brief This file implement the slicing visitors in the namespace STK.
 **/

/** @ingroup Arrays
 *  @defgroup SlicingVisitors Slicing Visitors
 * A slicing visitor is applied on each column/row of an Expression
 * or Array by using global functions in the STK domain space
 * <table class="example">
 * <tr><th>Example:</th><th>Output:</th></tr>
 * <tr><td>
 * \include tutoSlicingVisitors.cpp
 * </td>
 * <td>
 * \verbinclude tutoSlicingVisitors.out
 * </td></tr>
 * </table>
 **/


#ifndef STK_SLICEVISITORS_H
#define STK_SLICEVISITORS_H

#define STK_SLICEVISITORS(VISITOR, FUNC) \
template< class Derived> \
typename hidden::SliceVisitorSelector<Derived, hidden::VISITOR, Arrays::by_col_>::type_result \
FUNC(Derived const& A) \
{ return typename hidden::SliceVisitorSelector<Derived, hidden::VISITOR , Arrays::by_col_>::VisitorOp(A);} \
template< class Derived> \
typename hidden::SliceVisitorSelector<Derived, hidden::VISITOR, Arrays::by_col_>::type_result \
FUNC##ByCol(Derived const& A) \
{ return typename hidden::SliceVisitorSelector<Derived, hidden::VISITOR, Arrays::by_col_>::VisitorOp(A);} \
template< class Derived> \
typename hidden::SliceVisitorSelector<Derived, hidden::VISITOR, Arrays::by_row_>::type_result \
FUNC##ByRow(Derived const& A) \
{ return typename hidden::SliceVisitorSelector<Derived, hidden::VISITOR, Arrays::by_row_>::VisitorOp(A);}

#include "../allocators/STK_CAllocator.h"

namespace STK
{

// forward declaration
template<typename Derived, template<class> class Visitor>
class VisitorByCol;

namespace hidden
{

/** @ingroup hidden
 *  @brief Traits class for the unary operators
 */
template<typename Derived, template<class> class Visitor>
struct Traits< VisitorByCol <Derived, Visitor> >
{
  enum
  {
      structure_ = Arrays::point_,
      orient_    = Arrays::by_row_,
      sizeRows_  = 1,
      sizeCols_  = Derived::sizeCols_,
      storage_   = Arrays::dense_
  };
  typedef typename Derived::Type Type_;
  typedef Visitor<Type_> VisitorType;
  typedef typename VisitorType::Type Type;
  typedef typename VisitorType::ConstReturnType ConstReturnType;

//  typedef RowOperator< VisitorByCol <Derived, Visitor> > Row;
//  typedef ColOperator< VisitorByCol <Derived, Visitor> > Col;

  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

} // end namespace hidden

// forward declaration
template<typename Derived, template<class> class Visitor>
class VisitorByColBase;


/** @ingroup SlicingVisitors
 *  @brief  class allowing to apply a visitor on each columns of an expression.
 **/
template<typename Derived, template<class> class Visitor>
class VisitorByCol: public ExprBase< VisitorByCol<Derived, Visitor> >, public TRef<1>
{
  public:
    typedef typename hidden::Traits< VisitorByCol<Derived, Visitor> >::VisitorType VisitorType;
    typedef typename hidden::Traits< VisitorByCol<Derived, Visitor> >::Type Type;
    typedef typename hidden::Traits< VisitorByCol<Derived, Visitor> >::ConstReturnType ConstReturnType;
    typedef typename hidden::Traits< VisitorByCol<Derived, Visitor> >::Allocator Allocator;
    typedef VisitorByCol Result;
    enum
    {
      structure_ = hidden::Traits< VisitorByCol<Derived, Visitor> >::structure_,
      orient_    = hidden::Traits< VisitorByCol<Derived, Visitor> >::orient_,
      sizeRows_  = hidden::Traits< VisitorByCol<Derived, Visitor> >::sizeRows_,
      sizeCols_  = hidden::Traits< VisitorByCol<Derived, Visitor> >::sizeCols_,
      storage_   = hidden::Traits< VisitorByCol<Derived, Visitor> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** constructor */
    VisitorByCol( Derived const& lhs) : lhs_(lhs), result_(1,lhs_.sizeCols())
    {
      result_.shift(lhs_.beginCols());
      for (int j= lhs_.beginCols(); j < lhs_.endCols(); ++j)
      {
        VisitorType visit;
        result_.elt(j) = lhs_.col(j).visit(visit);
      }
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return result_.rows();}
    /** @return the columns range */
    inline ColRange const&colsImpl() const { return result_.cols();}

    /** @return the left hand side expression */
    inline Derived const& lhs() const { return lhs_; }
    /** @return the result */
    inline Allocator const& result() const { return result_; }

    /** access to the element (i,j) */
    inline Type const& elt2Impl(int i, int j) const { return result_.elt(i,j);}
    /** access to the element i */
    inline Type const& elt1Impl(int i) const { return result_.elt(i);}
    /** access to the element */
    inline Type const& elt0Impl() const { return result_.elt();}

  protected:
    Derived const& lhs_;

  private:
    Allocator result_;
};

// forward declaration
template<typename Derived, template<class> class Visitor>
class VisitorByRow;

namespace hidden
{

/** @ingroup hidden
 *  @brief Traits class for the VisitorByRow class
 */
template<typename Derived, template<class> class Visitor>
struct Traits< VisitorByRow <Derived, Visitor> >
{
  enum
  {
      structure_ = Arrays::vector_,
      orient_    = Arrays::by_col_,
      sizeRows_  = Derived::sizeRows_,
      sizeCols_  = 1,
      storage_   = Arrays::dense_
  };
  typedef typename Derived::Type Type_;
  typedef Visitor<Type_> VisitorType;
  typedef typename VisitorType::Type Type;
  typedef typename VisitorType::ConstReturnType ConstReturnType;

//  typedef RowOperator< VisitorByRow <Derived, Visitor> > Row;
//  typedef ColOperator< VisitorByRow <Derived, Visitor> > Col;

  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

} // end namespace hidden


/** @ingroup SlicingVisitors
 *  @brief  class allowing to apply a visitor on each columns of an expression.
 **/
template<typename Derived, template<class> class Visitor>
class VisitorByRow: public ExprBase< VisitorByRow<Derived, Visitor> >, public TRef<1>
{
  public:

    typedef typename hidden::Traits< VisitorByRow<Derived, Visitor> >::Type Type;
    typedef typename hidden::Traits< VisitorByRow<Derived, Visitor> >::VisitorType VisitorType;
    typedef typename hidden::Traits< VisitorByRow<Derived, Visitor> >::ConstReturnType ConstReturnType;
    typedef typename hidden::Traits< VisitorByRow<Derived, Visitor> >::Allocator Allocator;

    typedef VisitorByRow Result;
    enum
    {
      structure_ = hidden::Traits<VisitorByRow<Derived, Visitor> >::structure_,
      orient_    = hidden::Traits<VisitorByRow<Derived, Visitor> >::orient_,
      sizeRows_  = hidden::Traits<VisitorByRow<Derived, Visitor> >::sizeRows_,
      sizeCols_  = hidden::Traits<VisitorByRow<Derived, Visitor> >::sizeCols_,
      storage_   = hidden::Traits<VisitorByRow<Derived, Visitor> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** constructor */
    VisitorByRow( Derived const& lhs): lhs_(lhs), result_(lhs_.sizeRows(), 1)
    {
      result_.shift(lhs_.beginRows());
      for (int i= lhs_.beginRows(); i < lhs_.endRows(); ++i)
      {
        VisitorType visit;
        result_.elt(i) = lhs_.row(i).visit(visit);
      }
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return result_.rows();}
    /** @return the columns range */
    inline ColRange const&colsImpl() const { return result_.cols();}

    /** @return the left hand side expression */
    inline Derived const& lhs() const { return lhs_; }
    /** @return the result */
    inline Allocator const& result() const { return result_; }

    /** access to the element (i,j) */
    inline Type const& elt2Impl(int i, int j) const { return result_.elt(i,j);}
    /** access to the element i */
    inline Type const& elt1Impl(int i) const { return result_.elt(i);}
    /** access to the element */
    inline Type const& elt0Impl() const { return result_.elt();}

 protected:
    Derived const& lhs_;

  private:
    Allocator result_;
};

/** @ingroup SlicingVisitors
 *  @brief  class allowing to apply a visitor on a vector or a point.
 **/
template<typename Derived, template<class> class Visitor>
struct ApplyVisitor
{
  typedef typename Derived::Type Type_;
  typedef typename Visitor<Type_>::Type Type;
  typedef Type Result;

  /** constructor */
  inline ApplyVisitor( ExprBase<Derived> const& lhs) : lhs_(lhs.asDerived())
  {
    Visitor<Type_> visit;
    lhs_.visit(visit);
  }
  /** overload cast operator */
  inline operator Type () const { return result_;}
  protected:
    Derived const& lhs_;
  private:
     Type result_;
};

// Start utilites function
namespace hidden
{
/** @ingroup hidden
 *  Utility class that will select the type of operator to apply.
 *  The result can be either a number if the data are in a vector or a point,
 *  or a vector if the data are in a matrix
 **/
template<typename Derived, template<class> class Visitor, bool byCol>
struct SliceVisitorSelector;

/** @ingroup hidden
 *  Specialization if the visitor has to be applied by col
 **/
template<typename Derived, template<class> class Visitor>
struct SliceVisitorSelector< Derived, Visitor, true>
{
  typedef typename Derived::Type Type;
  enum
  {
    isVector_   =  (  Derived::structure_ == int(Arrays::vector_)
                   || Derived::structure_ == int(Arrays::number_)
                   )
  };
  typedef typename If<(bool)isVector_, ApplyVisitor<Derived, Visitor>, VisitorByCol<Derived, Visitor> >::Result VisitorOp;
  typedef typename VisitorOp::Result type_result;
};

/** @ingroup hidden
 *  Specialization if the visitor has to be applied by row
 **/
template<typename Derived, template<class> class Visitor>
struct SliceVisitorSelector< Derived, Visitor, false>
{
  typedef typename Derived::Type Type;
  enum
  {
    isPoint_   = (  Derived::structure_ == int(Arrays::point_)
                 || Derived::structure_ == int(Arrays::number_)
                 )
  };
  typedef typename If<(bool)isPoint_, ApplyVisitor<Derived, Visitor>, VisitorByRow<Derived, Visitor> >::Result VisitorOp;
  typedef typename VisitorOp::Result type_result;
};

} // hidden

/** @ingroup SlicingVisitors
 *  If A is a row-vector or a column-vector then the function will return the
 *  usual minimal value of the vector. If A is a two-dimensional array, the
 *  function will return a STK::VisitorByCol or a STK::VisitorByRow structure
 *  with the minimal values.
 *  @sa STK::max, STK::sum, STK::mean, STK::count, STK::all, STK::any.
 **/
STK_SLICEVISITORS(MinVisitor,min)
/** @ingroup SlicingVisitors
 *  If A is a row-vector or a column-vector then the function will return the
 *  usual maximal value of the vector. If A is a two-dimensional array, the
 *  function will return a STK::VisitorByCol or a STK::VisitorByRow structure
 *  with the maximal values.
 *  @sa STK::min, STK::sum, STK::mean, STK::count, STK::all, STK::any.
 **/
STK_SLICEVISITORS(MaxVisitor,max)
/** @ingroup SlicingVisitors
 *  If A is a row-vector or a column-vector then the function will return the
 *  usual sum value of the vector. If A is a two-dimensional array, the
 *  function will return a STK::VisitorByCol or a STK::VisitorByRow structure
 *  with the sums.
 *  @sa STK::min, STK::max, STK::mean, STK::count, STK::all, STK::any.
 **/
STK_SLICEVISITORS(SumVisitor,sum)
/** @ingroup SlicingVisitors
 *  If A is a row-vector or a column-vector then the function will return the
 *  usual mean value of the vector. If A is a two-dimensional array, the
 *  function will return a STK::VisitorByCol or a STK::VisitorByRow structure
 *  with the means.
 *  @sa STK::min, STK::max, STK::sum, STK::count, STK::all, STK::any,.
 **/
STK_SLICEVISITORS(MeanVisitor,mean)
/** @ingroup SlicingVisitors
 *  If A is a row-vector or a column-vector then the function will return the
 *  usual count value of the vector. If A is a two-dimensional array, the
 *  function will return a STK::VisitorByCol or a STK::VisitorByRow structure
 *  with the counts.
 *  @sa STK::min, STK::max, STK::sum, STK::mean, STK::all, STK::any.
 **/
STK_SLICEVISITORS(CountVisitor,count)
/** @ingroup SlicingVisitors
 *  If A is a row-vector or a column-vector then the function will return the
 *  usual all value of the vector. If A is a two-dimensional array, the
 *  function will return a STK::VisitorByCol or a STK::VisitorByRow structure
 *  with the alls values.
 *  @sa STK::min, STK::max, STK::sum, STK::mean, STK::count, STK::any.
 **/
STK_SLICEVISITORS(AllVisitor,all)
/** @ingroup SlicingVisitors
 *  If A is a row-vector or a column-vector then the function will return the
 *  usual visitor value of the vector. If A is a two-dimensional array, the
 *  function will return a STK::VisitorByCol or a STK::VisitorByRow structure
 *  with the any-s values.
 *  @sa STK::min, STK::max, STK::sum, STK::mean, STK::count, STK::all.
 **/
STK_SLICEVISITORS(AnyVisitor,any)

/** @ingroup SlicingVisitors
 *  If A is a row-vector or a column-vector then the function will return the
 *  usual minimal value of the vector. If A is a two-dimensional array, the
 *  function will return a STK::VisitorByCol or a STK::VisitorByRow structure
 *  with the minimal values.
 *  @sa STK::max, STK::sum, STK::mean, STK::count, STK::all, STK::any.
 **/
STK_SLICEVISITORS(MinSafeVisitor,minSafe)
/** @ingroup SlicingVisitors
 *  If A is a row-vector or a column-vector then the function will return the
 *  usual maximal value of the vector. If A is a two-dimensional array, the
 *  function will return a STK::VisitorByCol or a STK::VisitorByRow structure
 *  with the maximal values.
 *  @sa STK::min, STK::sum, STK::mean, STK::count, STK::all, STK::any.
 **/
STK_SLICEVISITORS(MaxSafeVisitor,maxSafe)
/** @ingroup SlicingVisitors
 *  If A is a row-vector or a column-vector then the function will return the
 *  usual mean value of the vector. If A is a two-dimensional array, the
 *  function will return a STK::VisitorByCol or a STK::VisitorByRow structure
 *  with the means.
 *  @sa STK::min, STK::max, STK::sum, STK::count, STK::all, STK::any,.
 **/
STK_SLICEVISITORS(MeanSafeVisitor,meanSafe)

}  // namespace STK

#undef STK_SLICEVISITORS

#endif /*STK_SLICEVISITORS_H*/
