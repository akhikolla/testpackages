/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

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
 * Project:  stkpp
 * created on: 28 juil. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file Rcpp_as.h
 *  @brief In this file we create the implement the non intrusive Rcpp::as Exporter class .
 **/

#ifndef RCPP_AS_H
#define RCPP_AS_H

namespace Rcpp
{

namespace traits
{
/* support for Rcpp::as */
template<typename Type>
class Exporter< STK::RMatrix<Type> >
{
  private:
    enum
    {
      Rtype_ = STK::hidden::RcppTraits<Type>::Rtype_
    };
    Rcpp::Matrix<Rtype_> mat;

  public:
    Exporter(SEXP x) : mat(x)
    {
      if (TYPEOF(x) != Rtype_) ::Rf_error("Wrong R type for mapped matrix");
    }
    STK::RMatrix<Type> get() {return STK::RMatrix<Type>(mat);}
} ;

} // namespace traits

} // namespace Rcpp

#endif // RCPP_AS_H
