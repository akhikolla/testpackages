/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2018  Serge Iovleff, Université Lille 1, Inria

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

/*    Project: stkpp::
 * created on: May 7, 2018
 *     Author: iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Functors.h
 *  @brief In this file we implement the distribution functors (pdf, lpdf, cdf, icdf) used by arrays
 **/



#ifndef STK_LAW_FUNCTORS_H
#define STK_LAW_FUNCTORS_H

#include "STK_Law_IUnivLaw.h"

namespace STK
{

namespace Law
{

/** @ingroup Law,Functors
 *  @brief Template functor computing  pdf values
 **/
template<typename Type>
struct PdfOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline PdfOp(IUnivLaw<Type> const& law): law_(law) {}
  inline PdfOp(PdfOp const& functor): law_(functor.law_) {}
  inline result_type operator()(param1_type const& a) const { return law_.pdf(a);}

  private:
    IUnivLaw<Type> const& law_;
};

/** @ingroup Law,Functors
 *  @brief Template functor computing  log-pdf values
 **/
template<typename Type>
struct LogPdfOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline LogPdfOp(IUnivLaw<Type> const& law): law_(law) {}
  inline LogPdfOp(LogPdfOp const& functor): law_(functor.law_) {}
  inline result_type operator()(param1_type const& a) const { return law_.lpdf(a);}

  private:
    IUnivLaw<Type> const& law_;
};

/** @ingroup Law,Functors
 *  @brief Template functor computing  cdf values
 **/
template<typename Type>
struct CdfOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline CdfOp(IUnivLaw<Type> const& law): law_(law) {}
  inline CdfOp(CdfOp const& functor): law_(functor.law_) {}
  inline result_type operator()(param1_type const& a) const { return law_.cdf(a);}

  private:
    IUnivLaw<Type> const& law_;
};

/** @ingroup Law,Functors
 *  @brief Template functor computing log-cdf values
 **/
template<typename Type>
struct LogCdfOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline LogCdfOp(IUnivLaw<Type> const& law): law_(law) {}
  inline LogCdfOp(LogCdfOp const& functor): law_(functor.law_) {}
  inline result_type operator()(param1_type const& a) const { return law_.lcdf(a);}

  private:
    IUnivLaw<Type> const& law_;
};

/** @ingroup Law,Functors
 *  @brief Template functor computing  complementary cdf values
 **/
template<typename Type>
struct CdfcOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline CdfcOp(IUnivLaw<Type> const& law): law_(law) {}
  inline CdfcOp(CdfcOp const& functor): law_(functor.law_) {}
  inline result_type operator()(param1_type const& a) const { return law_.cdfc(a);}

  private:
    IUnivLaw<Type> const& law_;
};

/** @ingroup Law,Functors
 *  @brief Template functor computing log-cdf complementary values
 **/
template<typename Type>
struct LogCdfcOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline LogCdfcOp(IUnivLaw<Type> const& law): law_(law) {}
  inline LogCdfcOp(LogCdfcOp const& functor): law_(functor.law_) {}
  inline result_type operator()(param1_type const& a) const { return law_.lcdfc(a);}

  private:
    IUnivLaw<Type> const& law_;
};

/** @ingroup Law,Functors
 *  @brief Template functor computing icdf values
 **/
template<typename Type>
struct IcdfOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline IcdfOp(IUnivLaw<Type> const& law): law_(law) {}
  inline IcdfOp(IcdfOp const& functor): law_(functor.law_) {}
  inline result_type operator()(param1_type const& a) const { return law_.icdf(a);}

  private:
    IUnivLaw<Type> const& law_;
};

} // namespace law

} // namespace SK

#endif /* STK_LAW_FUNCTORS_H */
