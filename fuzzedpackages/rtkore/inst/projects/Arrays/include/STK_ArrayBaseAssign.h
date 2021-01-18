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
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ArrayBaseAssign.h
 *  @brief In this file we implement the copy and assign methods
 **/


#ifndef STK_ARRAYBASEASSIGN_H
#define STK_ARRAYBASEASSIGN_H

#include "assign/STK_AssignSelector.h"

// this macro will be true if the assignation is correct and false otherwise
#define IS_VALID_ASSIGN(lhs,rhs) \
( (  ( lhs==Arrays::array2D_  || lhs==Arrays::square_) \
     && \
     ( rhs==Arrays::array2D_           || rhs==Arrays::square_  \
     || rhs==Arrays::diagonal_         || rhs==Arrays::number_ \
     || rhs==Arrays::lower_triangular_ || rhs==Arrays::upper_triangular_ \
     || rhs==Arrays::lower_symmetric_  || rhs==Arrays::upper_symmetric_ \
     || rhs==Arrays::symmetric_ \
     ) \
  )  \
  || \
  (  ( lhs==Arrays::array2D_) \
     && \
     (rhs==Arrays::vector_ || rhs==Arrays::point_|| rhs==Arrays::number_) \
  )  \
  || \
  ( lhs==Arrays::lower_triangular_ && rhs==Arrays::lower_triangular_) \
  || \
  ( lhs==Arrays::upper_triangular_ && rhs==Arrays::upper_triangular_) \
  || \
  ( ( lhs==Arrays::diagonal_ || lhs==Arrays::vector_ || lhs==Arrays::point_) \
    && \
    (rhs==Arrays::diagonal_ || rhs==Arrays::vector_ || rhs==Arrays::point_) \
  ) \
  || \
  ( lhs==Arrays::number_ && rhs==Arrays::number_) \
)

namespace STK
{

/* @brief assign rhs to this **/
template<class Derived>
template<class Rhs>
inline Derived& ArrayBase<Derived>::assign(ExprBase<Rhs> const& rhs)
{
  enum
  {
    rhs_structure_ = hidden::Traits<Rhs>::structure_
  , rhs_orient_    = hidden::Traits<Rhs>::orient_
  , rhs_sizeRows_  = hidden::Traits<Rhs>::sizeRows_
  , rhs_sizeCols_  = hidden::Traits<Rhs>::sizeCols_
  , is_valid_ = IS_VALID_ASSIGN((Arrays::Structure)structure_, (Arrays::Structure)rhs_structure_)
  };
  STK_STATIC_ASSERT(is_valid_,YOU_TRIED_TO_ASSIGN_A_NOT_COMPATIBLE_ARRAY);
  // check if assignment is possible
  if (structure_ == int(Arrays::square_) && rhs.cols() != rhs.rows() )
  { STKRUNTIME_ERROR_2ARG(ArrayBase::assign,Arrays::structureToString((Arrays::Structure)structure_),Arrays::structureToString((Arrays::Structure)rhs_structure_),is not permited);}
  // choose the correct way to resize this if necessary
  hidden::resizeSelector<Derived, Rhs, rhs_structure_>::run(this->asDerived(), rhs.asDerived());
  // choose the correct way to copy
  hidden::CopycatSelector<Derived, Rhs,  rhs_orient_>::run(this->asDerived(), rhs.asDerived());
  return this->asDerived();
}

/* @return the matrix or vector obtained by setting this constant*/
template<class Derived>
inline Derived& ArrayBase<Derived>::operator=( Type const& value) { return setValue(value);}
/* @return the matrix or vector obtained by evaluating this expression */
template<class Derived>
inline Derived& ArrayBase<Derived>::operator=( Derived const& rhs) { return assign(rhs);}
/* @return the matrix or vector obtained by evaluating this expression */
template<class Derived>
template<typename Rhs>
inline Derived& ArrayBase<Derived>::operator=( ExprBase<Rhs> const& rhs)
{ return assign(rhs.asDerived());}

//------------------------------------------------
// Rhs array
template<class Derived>
template<typename Rhs>
inline Derived&  ArrayBase<Derived>::operator+=( ExprBase<Rhs> const& rhs)
{
  enum { orient_     = hidden::Traits<Derived>::orient_
       , RStructure_ = hidden::Traits<Rhs>::structure_
       };
  typedef typename hidden::OperatorHelper<Derived, Rhs, Arrays::sumOp_>::Result Res;
  hidden::CopycatSelector<Derived, Res, orient_>::run(this->asDerived(), this->asDerived() + rhs.asDerived());
  return this->asDerived();
}
template<class Derived>
template<typename Rhs>
inline Derived&  ArrayBase<Derived>::operator-=( ExprBase<Rhs> const& rhs)
{
    enum { orient_     = hidden::Traits<Derived>::orient_
         , RStructure_ = hidden::Traits<Rhs>::structure_
         };
  typedef typename hidden::OperatorHelper<Derived, Rhs, Arrays::differenceOp_>::Result Res;
  hidden::CopycatSelector<Derived, Res, orient_>::run(this->asDerived(), this->asDerived() - rhs.asDerived());
  return this->asDerived();
}

template<class Derived>
template<typename Rhs>
inline Derived&  ArrayBase<Derived>::operator/=( ExprBase<Rhs> const& rhs)
{
    enum { orient_     = hidden::Traits<Derived>::orient_
         , RStructure_ = hidden::Traits<Rhs>::structure_
         };
  typedef typename hidden::OperatorHelper<Derived, Rhs, Arrays::divisionOp_>::Result Res;
  hidden::CopycatSelector<Derived, Res, orient_>::run(this->asDerived(), this->asDerived() / rhs.asDerived());
  return this->asDerived();
}
// TODO: implement direct call to CopyCat
template<class Derived>
template<typename Rhs>
inline Derived&  ArrayBase<Derived>::operator*=( ExprBase<Rhs> const& rhs)
{
    enum { orient_     = hidden::Traits<Derived>::orient_
         , RStructure_ = hidden::Traits<Rhs>::structure_
         };
  //typedef typename ProductProductType<Derived, Rhs>::ProductType Res;
//  typedef BinaryOperator< ProductOp<Type, typename hidden::Traits<Rhs>::Type>, Derived, Rhs> Res;
//  hidden::CopycatSelector<Derived, Res, orient_>::run(this->asDerived(), this->asDerived() * rhs.asDerived());
  this->asDerived() = this->asDerived() * rhs.asDerived();
  return this->asDerived();
}

template<class Derived>
template<typename Rhs>
inline Derived&  ArrayBase<Derived>::operator%=( ExprBase<Rhs> const& rhs)
{
    enum { orient_     = hidden::Traits<Derived>::orient_
         , RStructure_ = hidden::Traits<Rhs>::structure_
         };
  typedef typename hidden::OperatorHelper<Derived, Rhs, Arrays::moduloOp_>::Result Res;
  hidden::CopycatSelector<Derived, Res, orient_>::run(this->asDerived(), this->asDerived() % rhs.asDerived());
  return this->asDerived();
}

//------------------------------------------------
// rhs value
template<class Derived>
inline Derived& ArrayBase<Derived>::operator+=( Type const& value)
{
  enum { orient_ = hidden::Traits<Derived>::orient_};
  typedef UnaryOperator<SumWithOp<Type>, Derived> Rhs;
  hidden::CopycatSelector<Derived, Rhs, orient_>::run(this->asDerived(), this->asDerived() + value);
  return this->asDerived();
}
template<class Derived>
inline Derived& ArrayBase<Derived>::operator-=( Type const& value)
{
  enum { orient_ = hidden::Traits<Derived>::orient_};
  typedef UnaryOperator<SumWithOp<Type>, Derived> Rhs;
  hidden::CopycatSelector<Derived, Rhs, orient_>::run(this->asDerived(), this->asDerived() + (-value));
  return this->asDerived();
}
template<class Derived>
inline Derived& ArrayBase<Derived>::operator*=( Type const& value)
{
  enum { orient_ = hidden::Traits<Derived>::orient_};
  typedef UnaryOperator<ProductWithOp<Type>, Derived> Rhs;
  hidden::CopycatSelector<Derived, Rhs, orient_>::run(this->asDerived(), this->asDerived() * value);
  return this->asDerived();
}
template<class Derived>
inline Derived& ArrayBase<Derived>::operator/=( Type const& value)
{
  enum { orient_ = hidden::Traits<Derived>::orient_};
  typedef UnaryOperator<DivisionWithOp<Type>, Derived> Rhs;
  hidden::CopycatSelector<Derived, Rhs, orient_>::run(this->asDerived(), this->asDerived() / value);
  return this->asDerived();
}
template<class Derived>
inline Derived& ArrayBase<Derived>::operator%=( Type const& value)
{
  enum { orient_ = hidden::Traits<Derived>::orient_};
  typedef UnaryOperator<ModuloWithOp<Type>, Derived> Rhs;
  hidden::CopycatSelector<Derived, Rhs, orient_>::run(this->asDerived(), this->asDerived() % value);
  return this->asDerived();
}

/* overwrite @c this with @c rhs.
 *  @note this method does not take care of the possibility of overlapping
 *  @param rhs the array to copy
 **/
template<class Derived>
template<class Rhs>
inline Derived& ArrayBase<Derived>::copy( ExprBase<Rhs> const& rhs)
{
  // check
  if (this->sizeRows() != rhs.sizeRows())
  { STKRUNTIME_ERROR_2ARG(ArrayBase<Derived>::copy,this->sizeRows(), rhs.sizeRows(),sizeRows are not the sames);}
  if (this->sizeCols() != rhs.sizeCols())
  { STKRUNTIME_ERROR_2ARG(ArrayBase<Derived>::copy,this->sizeCols(), rhs.sizeCols(),sizeCols are not the sames);}
  // copy
  for ( int jRhs=rhs.beginCols(), jLhs=this->beginCols(); jRhs<rhs.endCols(); jLhs++, jRhs++)
    for ( int iRhs=rhs.beginRows(), iLhs=this->beginRows(); iRhs<rhs.endRows(); iLhs++, iRhs++)
  { this->elt(iLhs, jLhs) = rhs(iRhs, jRhs);}
  // return this
  return this->asDerived();
}

} // namespace STK

#undef IS_VALID_ASSIGN

#endif /* STK_ARRAYBASEASSIGN_H */
