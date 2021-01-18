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

/** @file STK_ExprBaseFunctors.h
 *  @brief In this file we define the way to apply functors to Expressions.
 **/

#ifndef STK_EXPRBASEFUNCTORS_H
#define STK_EXPRBASEFUNCTORS_H

#include "STK_ExprBase.h"

namespace STK
{

// forward declaration
template<typename Derived, typename Functor> struct ApplyFunctorByCol;
template<typename Derived, typename Functor> struct ApplyFunctorByRow;
template<typename Derived, typename Functor> struct ApplyFunctor;
template<typename Derived, typename Functor> struct ApplyWeightedFunctorByCol;
template<typename Derived, typename Functor> struct ApplyWeightedFunctorByRow;
template<typename Derived, typename Functor> struct ApplyWeightedFunctor;

namespace hidden
{
/** @ingroup hidden, Array
 *  @brief Utility class allowing to get the correct Structure, Type, allocator,...
 *  of the result when a functor is applied to an array.
 *
 *  @sa MinOp, MaxOp, MeanOp,MeanSafeOp
 **/
template<typename Derived, template<class> class Functor, bool isVector_>
struct FunctorSelector;


/** @ingroup hidden
 *  specialization for general arrays: apply functor to rows or columns
 *  of the @c Derived array
 **/
template<class Derived, template<class> class Functor>
struct FunctorSelector<Derived, Functor, false>
{
  //typedef typename Derived::Col Col;
  typedef typename hidden::Traits<Derived>::Col Col;
  typedef Functor<Col> ByColFunctor;
  typedef ApplyFunctorByCol<Derived, ByColFunctor> ColOp;
  typedef ApplyWeightedFunctorByCol<Derived, ByColFunctor> ColWeightedOp;

  //typedef typename Derived::Row Row;
  typedef typename hidden::Traits<Derived >::Row Row;
  typedef Functor<Row> ByRowFunctor;
  typedef ApplyFunctorByRow<Derived, ByRowFunctor> RowOp;
  typedef ApplyWeightedFunctorByRow<Derived, ByRowFunctor> RowWeightedOp;
};

/** @ingroup hidden
 *  specialization for vectors and points. In this case the functor apply
 *  directly to the @c Derived class.
 **/
template<class Derived, template<class> class Functor>
struct FunctorSelector<Derived, Functor, true>
{
  typedef Functor<Derived> ByColFunctor;
  typedef ApplyFunctor<Derived, ByColFunctor> ColOp;
  typedef ApplyWeightedFunctor<Derived, ByColFunctor> ColWeightedOp;

  typedef Functor<Derived> ByRowFunctor;
  typedef ApplyFunctor<Derived, ByRowFunctor> RowOp;
  typedef ApplyWeightedFunctor<Derived, ByRowFunctor> RowWeightedOp;
};

/** @ingroup hidden
 *  Utility class that will select the type of operator to apply.
 *  The result can be either a number if the data are in a vector or a point,
 *  or a vector if the data are in a matrix
 **/
template<typename Derived, template<class> class Functor>
struct FunctorTraits
{
  enum
  {
    isVector_   =  (  Derived::structure_ == int(Arrays::vector_)
                   || Derived::structure_ == int(Arrays::point_)
                   || Derived::structure_ == int(Arrays::number_)
                   )
  };

  typedef typename FunctorSelector<Derived, Functor, (bool)isVector_>::ColOp ColOp;
  typedef typename FunctorSelector<Derived, Functor, (bool)isVector_>::ColWeightedOp ColWeightedOp;
  typedef typename FunctorSelector<Derived, Functor, (bool)isVector_>::RowOp RowOp;
  typedef typename FunctorSelector<Derived, Functor, (bool)isVector_>::RowWeightedOp RowWeightedOp;

  typedef typename ColOp::Row Row;
  typedef typename RowOp::Col Col;
};

} // namespace hidden

/** @ingroup Arrays
 *  @brief  class allowing to apply the Functor Functor on each columns of an expression.
 **/
template<typename Derived, typename Functor>
struct ApplyFunctorByCol
{
  typedef typename Derived::Type Type;
  enum
  {
    sizeRows_ = Derived::sizeRows_,
    sizeCols_ = Derived::sizeCols_
  };
  typedef CArrayPoint<Type, sizeCols_> Row;

  /** constructor */
  inline ApplyFunctorByCol( ExprBase<Derived> const& lhs)
                          : lhs_(lhs.asDerived())
                          , res_(lhs_.cols())
                          {}
  /** @return the applied functor by column */
  Row operator()()
  {
    for (int j= lhs_.beginCols(); j < lhs_.endCols(); ++j)
    { res_[j] = Functor(lhs_.col(j))();}
    return res_;
  }
  /** @return the applied functor by column */
  Row operator()(bool option)
  {
    for (int j= lhs_.beginCols(); j < lhs_.endCols(); ++j)
    { res_[j] = Functor(lhs_.col(j))( option);}
    return res_;
  }
  /** @return the applied functor by column */
  template<class OtherRow>
  Row operator()(OtherRow const& value, bool option)
  {
#ifdef STK_BOUNDS_CHECK
  if (lhs_.cols() != value.range())
    STKRUNTIME_ERROR_NO_ARG(ApplyFunctorByCol::operator(value,option),lhs_.cols()!=value.range());
#endif
    for (int j= lhs_.beginCols(); j < lhs_.endCols(); ++j)
    { res_[j] = Functor(lhs_.col(j))(value[j], option);}
    return res_;
  }
  protected:
    Derived const& lhs_;
    Row res_;
};

/** @ingroup Arrays
 *  @brief  class allowing to apply the weighted Functor @c Functor on each
 *  columns of an expression.
 **/
template<typename Derived, typename Functor>
struct ApplyWeightedFunctorByCol
{
  typedef typename Derived::Type Type;
  enum
  {
    sizeRows_ = Derived::sizeRows_,
    sizeCols_ = Derived::sizeCols_
  };
  typedef CArrayPoint<Type, sizeCols_> Row;

  /** constructor */
  inline ApplyWeightedFunctorByCol( ExprBase<Derived> const& lhs)
                                  : lhs_(lhs.asDerived())
                                  , res_(lhs_.cols())
                                  {}
  template<class Weights>
  Row operator()(ExprBase<Weights> const& w)
  {
    for (int j= lhs_.beginCols(); j < lhs_.endCols(); ++j)
    { res_[j] = Functor(lhs_.col(j))(w);}
    return res_;
  }
  template<class Weights>
  Row operator()(ExprBase<Weights> const& w, bool option)
  {
    for (int j= lhs_.beginCols(); j < lhs_.endCols(); ++j)
    { res_[j] = Functor(lhs_.col(j))(w, option);}
    return res_;
  }
  template<class Weights>
  Row operator()(ExprBase<Weights> const& w, Type const& value)
  {
    for (int j= lhs_.beginCols(); j < lhs_.endCols(); ++j)
    { res_[j] = Functor(lhs_.col(j))(w, value);}
    return res_;
  }
  template< class Weights, typename OtherRow>
  Row operator()(ExprBase<Weights> const& w, ExprBase<OtherRow> const& value, bool option)
  {
#ifdef STK_BOUNDS_CHECK
  if (lhs_.cols() != value.range())
    STKRUNTIME_ERROR_NO_ARG(ApplyWeightedFunctorByCol::operator(w,value,option),lhs_.cols()!=value.range());
#endif
    for (int j= lhs_.beginCols(); j < lhs_.endCols(); ++j)
    { res_[j] = Functor(lhs_.col(j))(w, value[j], option);}
    return res_;
  }
  protected:
    Derived const& lhs_;
    Row res_;
};

/** @ingroup Arrays
 *  @brief class allowing to apply the Functor Functor on each rows of an expression.
 **/
template<typename Derived, typename Functor>
struct ApplyFunctorByRow
{
  typedef typename Derived::Type Type;
  enum
  {
    sizeRows_ = Derived::sizeRows_,
    sizeCols_ = Derived::sizeCols_
  };
  typedef CArrayVector<Type, sizeRows_> Col;

  /** constructor */
  inline ApplyFunctorByRow( ExprBase<Derived> const& lhs)
                          : lhs_(lhs.asDerived())
                          , res_(lhs_.rows())
                          {}
  Col operator()()
  {
    for (int i= lhs_.beginRows(); i < lhs_.endRows(); ++i)
    { res_[i] = Functor(lhs_.row(i))();}
    return res_;
  }
  Col operator()(bool option)
  {
    for (int j= lhs_.beginRows(); j < lhs_.endRows(); ++j)
    { res_[j] = Functor(lhs_.row(j))( option);}
    return res_;
  }
  template<class OtherCol>
  Col operator()(OtherCol const& value, bool option)
  {
#ifdef STK_BOUNDS_CHECK
  if (lhs_.rows() != value.range()) STKRUNTIME_ERROR_NO_ARG(ApplyFunctorByRow::operator(value,option),lhs_.rows()!=value.rows());
#endif
    for (int j= lhs_.beginRows(); j < lhs_.endRows(); ++j)
    { res_[j] = Functor(lhs_.row(j))(value[j], option);}
    return res_;
  }
  protected:
    Derived const& lhs_;
    Col res_;
};

/** @ingroup Arrays
 *  @brief class allowing to apply the Functor Functor on each rows of an expression.
 **/
template<typename Derived, typename Functor>
struct ApplyWeightedFunctorByRow
{
  typedef typename Derived::Type Type;
  enum
  {
    sizeRows_ = Derived::sizeRows_,
    sizeCols_ = Derived::sizeCols_
  };
  typedef CArrayVector<Type, sizeRows_> Col;

  /** constructor */
  inline ApplyWeightedFunctorByRow( ExprBase<Derived> const& lhs)
                                  : lhs_(lhs.asDerived())
                                  , res_(lhs_.rows())
                                  {}
  template<typename Weights>
  Col operator()( ExprBase<Weights> const& w)
  {
    for (int i= lhs_.beginRows(); i < lhs_.endRows(); ++i)
    { res_[i] = Functor(lhs_.row(i))(w);}
    return res_;
  }
  template<typename Weights>
  Col operator()(ExprBase<Weights> const& w, bool option)
  {
    for (int j= lhs_.beginRows(); j < lhs_.endRows(); ++j)
    { res_[j] = Functor(lhs_.row(j))(w, option);}
    return res_;
  }
  template<typename Weights, typename OtherCol>
  Col operator()(ExprBase<Weights> const& w, ExprBase<OtherCol> const& value, bool option)
  {
#ifdef STK_BOUNDS_CHECK
  if (lhs_.rows() != value.range()) STKRUNTIME_ERROR_NO_ARG(ApplyFunctorByRow::operator(w,value,option),lhs_.rows()!=value.rows());
#endif
    for (int j= lhs_.beginRows(); j < lhs_.endRows(); ++j)
    { res_[j] = Functor(lhs_.row(j))(w, value[j], option);}
    return res_;
  }
  protected:
    Derived const& lhs_;
    Col res_;
};

/** @ingroup Arrays
 *  @brief class allowing applying the functor @c Functor on a vector or row-vector
 **/
template<typename Derived, typename Functor>
struct ApplyFunctor
{
  typedef typename Derived::Type Type;
  typedef Type Col;
  typedef Type Row;
  /// constructor
  inline ApplyFunctor( ExprBase<Derived> const& lhs) : lhs_(lhs.asDerived())
  { STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);}
  /** apply without argument*/
  inline Type operator()() { return Functor(lhs_)();}
  /** apply with an option argument*/
  inline Type operator()(bool option) { return Functor(lhs_)(option);}
  /** apply with a value and an option argument*/
  inline Type operator()(Type const& value, bool option)
  { return Functor(lhs_)(value, option);}

  protected:
    Derived const& lhs_;
};

/** @ingroup Arrays
 *  @brief class allowing applying the weighted functor @c Functor on a vector or row-vector
 **/
template<typename Derived, typename Functor>
struct ApplyWeightedFunctor
{
  typedef typename Derived::Type Type;
  typedef Type Col;
  typedef Type Row;
  /// constructor
  inline ApplyWeightedFunctor( ExprBase<Derived> const& lhs): lhs_(lhs.asDerived())
  { STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);}
  /** apply with weights*/
  template<typename Weights>
  inline Type operator()(ExprBase<Weights> const& w) { return Functor(lhs_)(w);}
  /** apply with weight and an option argument*/
  template<typename Weights>
  inline Type operator()(ExprBase<Weights> const& w, bool option)
  { return Functor(lhs_)(w, option);}
  /** apply with weight, a value and an option argument*/
  template<typename Weights>
  inline Type operator()(ExprBase<Weights> const& w, Type const& value, bool option)
  { return Functor(lhs_)(w, value, option);}

  protected:
    Derived const& lhs_;
};


} // namespace STK

#endif /* STK_EXPRBASEVISITOR_H */
