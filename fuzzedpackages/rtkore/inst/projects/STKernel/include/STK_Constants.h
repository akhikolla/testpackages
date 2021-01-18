/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff

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
 * Project:  stkpp::stkernel
 * created on: 14 sept. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Constants.h
 *  @brief In this file we define the main constant which will be used through
 *  the STK++ project.
 **/

#ifndef STK_CONSTANTS_H
#define STK_CONSTANTS_H

#include <climits>

#ifndef STKBASEARRAYS
#define STKBASEARRAYS 0 // default is 0 based arrays
#endif

#ifndef MAXUNROLLVALUE
#define MAXUNROLLVALUE 20 // unroll operations on fixed size containers
#endif

#ifndef MAXUNROLLSLICEVALUE
#define MAXUNROLLSLICEVALUE 5 // unroll operations rows or columns on fixed size 2D containers
#endif

#ifndef MAXFIXEDSIZEVALUE
#define MAXFIXEDSIZEVALUE 1024 // maximal fixed size that can be used automatically
#endif


namespace STK
{
/** @ingroup STKernel
 *  @brief base index of the containers created in STK++.
 *  This value means that the default range for a vector or the rows/columns of
 *  a matrix is the value given by this constant. **/
const int baseIdx = STKBASEARRAYS;

/** @ingroup STKernel
 *  @brief maximal size of fixed size containers
 *  This value is used when fixed size containers are automatically build in internal computation
 **/
const int maxFixedSize = MAXFIXEDSIZEVALUE;

/** @ingroup STKernel
 *  This value means that an integer is not known at compile-time, and that
 *  instead the value is stored in some runtime variable. This is the same value
 *  that the value used for representing NA Integers when Integer is int.
 **/
const int UnknownSize = INT_MAX;
/** @ingroup STKernel
 *  Same as floor(sqrt(INT_MAX+1))
 **/
const int SqrtUnknownSize = (1 << (sizeof(int) * (CHAR_BIT/2)));

/** @ingroup STKernel
 *  When don't unroll loops on fixed size containers if fixed size is greater than this value.
 **/
const int MaxUnroll = MAXUNROLLVALUE;

/** @ingroup STKernel
 * This value means that when we unroll loops we go until MaxUnrollSlice
 **/
const int MaxUnrollSlice = MAXUNROLLSLICEVALUE;

} // namespace STK


#endif /* STK_CONSTANTS_H */
