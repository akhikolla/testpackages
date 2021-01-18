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
 * Project:  Algebra
 * Purpose:  Define the transpose method
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_transpose.h
 *  @brief In this file we define the transpose method.
 *
 **/

#ifndef STK_TRANSPOSE_H
#define STK_TRANSPOSE_H

// Container classes
#include <Arrays/include/STK_IArrayBase.h>

namespace STK
{

/** @ingroup Algebra
 *  @brief The transpose method allow to transpose an array in place
 *  @param A the container to transpose
 **/
template < class TContainer2D>
TContainer2D& transpose( ArrayBase<TContainer2D>& A)
{
  if(A.rows()!=A.cols())
  {
    TContainer2D Aux = A.transpose();
    A = Aux;
  }
  else // transpose in place
  {
    for (int i= A.beginRows(); i< A.endRows(); i++)
      for (int j= i+1; j< A.endCols(); j++) { std::swap(A(i,j),A(j,i));}
  }
  return A.asDerived();
}

} // namespace STK

#endif // STK_TRANSPOSE_H

