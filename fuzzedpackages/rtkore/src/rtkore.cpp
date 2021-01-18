/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015  Serge Iovleff

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
 * Project:  rtkore
 * created on: 1 mars 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file rtkore.cpp
 *  @brief In this file .
 **/


// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppEigen.cpp: Rcpp/Eigen glue
//
// Copyright (C)       2011 Douglas Bates, Dirk Eddelbuettel and Romain Francois
//
// This file is part of RcppEigen.
//
// RcppEigen is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppEigen is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppEigen.  If not, see <http://www.gnu.org/licenses/>.

#include "RTKpp.h"

using namespace STK;

RcppExport  SEXP stk_version(SEXP single)
{
  using Rcpp::_;
  using Rcpp::IntegerVector;

  BEGIN_RCPP;
  bool csingle = Rcpp::as<bool>(single) ;
  if( csingle)
  {
    return Rcpp::wrap(10000 * STK_WORLD_VERSION
                     +  100 * STK_MAJOR_VERSION
                     +        STK_MINOR_VERSION);
  }
  // return a list with the number versions
  return IntegerVector::create(_["major"] = STK_WORLD_VERSION,
                               _["minor"] = STK_MAJOR_VERSION,
                               _["patch"] = STK_MINOR_VERSION);
  END_RCPP;
}







