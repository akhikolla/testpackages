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

/** @file createKernel.cpp
 *  @brief In this file we create a kernel and optionally compute the gram matrix
 **/


#include "../inst/projects/MixAll/MixAll_Util.h"

/* Compute the gram matrix and store it in the "gram" slot of the component
 *  @param s4_component  the component with the data storing the gram matrix
 *  @param r_kernelName name of the kernel to use
 *  @param r_kernelParameters vector of parameters
 **/
STK::Kernel::IKernel* createKernel( Rcpp::S4 s4_component
                                  , STK::String const& kernelName
                                  , Rcpp::DoubleVector const& r_kernelParameters
                                  , bool computeGramMatrix
                                  )
{
  STK::Real param1, param2;
  switch (r_kernelParameters.length())
  {
    case 0:
      param1 = 1.; param2 = 0.;
      break;
    case 1:
      param1 = r_kernelParameters[0]; param2 = 0.;
      break;
    default:
      param1 = r_kernelParameters[0]; param2 =  r_kernelParameters[1];
      break;
  }
  STK::Kernel::IKernel* p_kernel =0;

  STK::RMatrix<double> datad;
  STK::RMatrix<int> datai;
  // build gram matrix and overwrite solt data with it
  switch (STK::Kernel::stringToKernelType(kernelName))
  {
    case STK::Kernel::laplace_:
      datad = s4_component.slot("data");
      p_kernel = new STK::Kernel::Laplace<STK::RMatrix<double> >(datad, param1);
      break;
    case STK::Kernel::gaussian_:
      datad = s4_component.slot("data");
      p_kernel = new STK::Kernel::Gaussian<STK::RMatrix<double> >(datad, param1);
      break;
    case STK::Kernel::linear_:
      datad = s4_component.slot("data");
      p_kernel = new STK::Kernel::Linear<STK::RMatrix<double> >(datad);
      break;
    case STK::Kernel::polynomial_:
      datad = s4_component.slot("data");
      p_kernel = new STK::Kernel::Polynomial<STK::RMatrix<double> >(datad, param1, param2);
      break;
    case STK::Kernel::rationalQuadratic_:
      datad = s4_component.slot("data");
      p_kernel = new STK::Kernel::RationalQuadratic<STK::RMatrix<double> >(datad, param1);
      break;
    case STK::Kernel::hamming_:
      datai = s4_component.slot("data");
      p_kernel = new STK::Kernel::Hamming<STK::RMatrix<int> >(datai, param1);
      break;
    default:
      break;
  }
  if (!p_kernel) return 0;
  // compute Gram matrix if needed
  if (computeGramMatrix)
  {
    if (!p_kernel->run()) { delete p_kernel; return 0;}
    s4_component.slot("gram") = STK::wrap(p_kernel->gram());
  }
  // return result
  return p_kernel;
}

