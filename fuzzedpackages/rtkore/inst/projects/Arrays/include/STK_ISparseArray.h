/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2017  Serge Iovleff, Université Lille 1, Inria

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
 * created on: 10 mars 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ISparseArray.h
 *  @brief In this file we define the base class for Sparse Matrices.
 **/

#ifndef STK_ISPARSEMATRIX_H
#define STK_ISPARSEMATRIX_H

#include <vector>
#include <string>

namespace STK
{

/** @class ISparseMatrix
  * @ingroup Arrays
  *
  * @brief Interface class for SparseArray, SparseVector, SparsePoint.
  *
  * This class is the base class that is inherited by all objects (matrix, vector,
  * point) which are not expressions and stored as SparseArrays. The common API
  * for these objects is contained in this class.
  *
  * @tparam Derived is the derived type, e.g., a matrix type.
  **/
template<class Derived>
class STK_ISparseArray: public ArrayBase<Derived>
{
  public:
    typedef ArrayBase<Derived> Base;

    typedef typename hidden::Traits<Derived>::Type Type;
    //typedef typename hidden::Traits<Derived>::Row  Row;
    typedef typename hidden::Traits<Derived>::Col  Col;
    //typedef typename hidden::Traits<Derived>::SubRow SubRow;
    //typedef typename hidden::Traits<Derived>::SubCol SubCol;
    //typedef typename hidden::Traits<Derived>::SubVector SubVector;
    typedef typename hidden::Traits<Derived>::SubArray SubArray;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    typedef hidden::CheckShift<Derived, structure_> CheckShift;
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    STK_ISparseArray ()
    {}

    ~STK_ISparseArray() {}

};

}

#endif /* STK_ISPARSEMATRIX_H */
