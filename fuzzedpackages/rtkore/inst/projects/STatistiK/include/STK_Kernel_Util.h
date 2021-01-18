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
 * Project:  stkpp::STatitiK
 * created on: 5 avr. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Kernel_Util.h
 *  @brief In this file we define the enum and utilites method for the kernels.
 **/


#ifndef STK_KERNEL_UTIL_H
#define STK_KERNEL_UTIL_H

#include <STKernel/include/STK_String.h>

namespace STK
{

namespace Kernel
{
/** @ingroup Kernel
 *  @brief kernel types.
 **/
enum kernelType
{
  laplace_,
  exponential_, // deprecated
  gaussian_,
  linear_,
  polynomial_,
  rationalQuadratic_,
  hamming_,
  unknown_kernel_
};

/** @ingroup Kernel
 *  Convert a String to a kernelType. The recognized strings are
 * <table >
 * <tr> <th> Kernel             </th> </tr>
 * <tr> <td> "Laplace"   </td></tr>
 * <tr> <td> "Gaussian"    </td></tr>
 * <tr> <td> "Linear"    </td></tr>
 * <tr> <td> "Polynomial"     </td></tr>
 * <tr> <td> "RationalQuadratic"    </td></tr>
 * <tr> <td> "Hamming"    </td></tr>
 * </table>
 *  @param type the String we want to convert
 *  @return the kernrlType represented by the String @c type. If the string
 *  does not match any known name, the @c unknown_kernel_ type is returned.
 **/
inline kernelType stringToKernelType( std::string const& type)
{
  if (toUpperString(type) == toUpperString(_T("Laplace"))) return laplace_;
  if (toUpperString(type) == toUpperString(_T("Exponential"))) return exponential_;
  if (toUpperString(type) == toUpperString(_T("Gaussian"))) return gaussian_;
  if (toUpperString(type) == toUpperString(_T("Linear"))) return linear_;
  if (toUpperString(type) == toUpperString(_T("Polynomial"))) return polynomial_;
  if (toUpperString(type) == toUpperString(_T("RationalQuadratic"))) return rationalQuadratic_;
  if (toUpperString(type) == toUpperString(_T("Hamming"))) return hamming_;
  return unknown_kernel_;
}

/** @ingroup Kernel
 *  convert a kernelType to a String.
 *  @param type the type of kernelType we want to convert
 *  @return the string associated to this type.
 **/
inline String kernelTypeToString( kernelType const& type)
{
  if (type == exponential_) return String(_T("Exponential")); // deprecated
  if (type == laplace_) return String(_T("Laplace"));
  if (type == gaussian_) return String(_T("Gaussian"));
  if (type == linear_) return String(_T("Linear"));
  if (type == polynomial_) return String(_T("Polynomial"));
  if (type == rationalQuadratic_) return String(_T("RationalQuadratic"));
  if (type == hamming_) return String(_T("Hamming"));
  return String(_T("unknown"));
}

} // namespace Kernel

} // namespace STK

#endif /* STK_KERNEL_UTIL_H */
