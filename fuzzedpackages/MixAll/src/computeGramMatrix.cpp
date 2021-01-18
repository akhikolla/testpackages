/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 28 juil. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file kmmMixAll.cpp
 *  @brief In this file we launch the computation for estimating a kernel mixture model.
 **/


#include "RTKpp.h"
#include "MixAll.h"
#include <MixAll/MixAll_Util.h>

/* Compute the Gram matrix and overwrite the data with the result.
 *  @param component a KmmModelComponent S4 class
 *  @param kernelName a string with the name of the kernel to use
 *  @param kernelParameters a vector with the optional parameters
 */
extern "C" SEXP computeGramMatrix( SEXP component, SEXP kernelName, SEXP kernelParameters)
{
  BEGIN_RCPP

  Rcpp::S4 s4_component(component);
  Rcpp::CharacterVector r_kernelName(kernelName);
  Rcpp::DoubleVector r_kernelParameters(kernelParameters);
  return(Rcpp::wrap(computeKernel(s4_component, r_kernelName, r_kernelParameters)));

  END_RCPP
}

