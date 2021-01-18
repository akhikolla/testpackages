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
 * Project:  MixAll::
 * created on: 23 juil. 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file MixAll_Util.h
 *  @brief In this file .
 **/


#ifndef MIXALL_UTIL_H_
#define MIXALL_UTIL_H_

#include "RTKpp.h"

namespace STK
{
// forward declaration
class IMixtureCriterion;
}
/** Compute the Gram matrix and overwrite the data with the result.
 *  @param s4_component a KmmComponent S4 class
 *  @param r_kernelName a string with the name of the kernel to use
 *  @param r_kernelParameters a vector with the optional parameters
 */
bool computeKernel( Rcpp::S4 s4_component
                  , Rcpp::CharacterVector const& r_kernelName
                  , Rcpp::DoubleVector const& r_kernelParameters
                  );

/** Create the kernel needed by the KMM model and compute the Gram matrix if needed
 *  @param s4_component  the component with the data storing the gram matrix
 *  @param r_kernelName name of the kernel to use
 *  @param r_kernelParameters vector of parameters
 *  @param computeGramMatrix @c true if the Gram matrix has to be computed, @c false
 *  otherwise
 **/
STK::Kernel::IKernel* createKernel( Rcpp::S4 s4_component
                                  , STK::String const& r_kernelName
                                  , Rcpp::DoubleVector const& r_kernelParameters
                                  , bool computeGramMatrix
                                  );


#endif /* MIXALL_UTIL_H_ */
