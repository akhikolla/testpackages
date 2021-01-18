
/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2015  <MODAL team @INRIA,Lille & U.M.R. C.N.R.S. 6599 Heudiasyc, UTC>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : parmeet.bhatia@inria.fr , bhatia.parmeet@gmail.com
*/


/** @file typedef.h
 *  @brief This file define all the typedefs used in cocluster project.
 **/


#ifndef TYPEDEF_H_
#define TYPEDEF_H_

#include <limits>
#ifdef RPACKAGE
#include <RTKpp.h>
#else
#include <STKpp.h>
#endif

//class InputParameters;
class ICoClustModel;
/*
 * Macro definition for numeric limits
 */
#define RealMax std::numeric_limits<STK::Real>::max()
#define RealMin std::numeric_limits<STK::Real>::min()

/*
 * Typedefs for Matrix and vector containers
 */

//Matrix containers
typedef STK::CArrayXX MatrixReal;
typedef STK::CArrayXXi MatrixInt;
typedef STK::CArray<bool, STK::UnknownSize, STK::UnknownSize> MatrixBinary;

//Vector Containers
typedef STK::CVectorXd VectorReal;
typedef STK::CPointXd PointReal;
typedef STK::CVectorXi VectorInt;
typedef STK::CArrayVector<bool, STK::UnknownSize> VectorBinary;

typedef STK::RVector<int> RVectorInt;

//2D array containers
typedef STK::ArrayXX Array2DReal;

/**
 * Member Function pointers
 */

typedef void(ICoClustModel::*StopCriteria_poiter)();

#endif /* TYPEDEF_H_ */
