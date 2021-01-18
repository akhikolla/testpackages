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

/** @file clusterMixture.cpp
 *  @brief In this file we launch the computation for estimating a mixture model.
 **/


#include <RTKpp.h>
#include <MixAll.h>

/* Compute the Gram matrix and overwrite the data with the result.
 *  @param s4_component a ClusterKernelComponent S4 class
 *  @param r_kernelName a string with the name of the kernel to use
 *  @param r_kernelParameters a vector with the optional parameters
 */
bool computeKernel( Rcpp::S4 s4_component
                  , Rcpp::CharacterVector const& r_kernelName
                  , Rcpp::DoubleVector const& r_kernelParameters
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
  STK::String kernelName = Rcpp::as<STK::String>(r_kernelName[0]);
  STK::Kernel::IKernelBase<STK::RMatrix<double> >* p_kerneld =0;
  STK::Kernel::IKernelBase<STK::RMatrix<int> >* p_kerneli =0;
  STK::RMatrix<double> datad;
  STK::RMatrix<int> datai;
  // build gram matrix and overwrite solt data with it
  switch (STK::Kernel::stringToKernelType(kernelName))
  {
    case STK::Kernel::exponential_:
      datad = s4_component.slot("rawData");
      p_kerneld = new STK::Kernel::Exponential<STK::RMatrix<double> >(datad, param1);
      if (!p_kerneld->run()) { delete p_kerneld; return Rcpp::wrap(false);}
      s4_component.slot("data") = STK::wrap(p_kerneld->gram());
      break;
    case STK::Kernel::gaussian_:
      datad = s4_component.slot("rawData");
      p_kerneld = new STK::Kernel::Gaussian<STK::RMatrix<double> >(datad, param1);
      if (!p_kerneld->run()) { delete p_kerneld; return Rcpp::wrap(false);}
      s4_component.slot("data") = STK::wrap(p_kerneld->gram());
      break;
    case STK::Kernel::linear_:
      datad = s4_component.slot("rawData");
      p_kerneld = new STK::Kernel::Linear<STK::RMatrix<double> >(datad);
      if (!p_kerneld->run()) { delete p_kerneld; return Rcpp::wrap(false);}
      s4_component.slot("data") = STK::wrap(p_kerneld->gram());
      break;
    case STK::Kernel::polynomial_:
      datad = s4_component.slot("rawData");
      p_kerneld = new STK::Kernel::Polynomial<STK::RMatrix<double> >(datad, param1, param2);
      if (!p_kerneld->run()) { delete p_kerneld; return Rcpp::wrap(false);}
      s4_component.slot("data") = STK::wrap(p_kerneld->gram());
      break;
    case STK::Kernel::rationalQuadratic_:
      datad = s4_component.slot("rawData");
      p_kerneld = new STK::Kernel::RationalQuadratic<STK::RMatrix<double> >(datad, param1);
      if (!p_kerneld->run()) { delete p_kerneld; return Rcpp::wrap(false);}
      s4_component.slot("data") = STK::wrap(p_kerneld->gram());
      break;
    case STK::Kernel::hamming_:
      datai = s4_component.slot("rawData");
      p_kerneli = new STK::Kernel::Hamming<STK::RMatrix<int> >(datai, param1);
      if (!p_kerneli->run()) { delete p_kerneli; return Rcpp::wrap(false);}
      s4_component.slot("data") = STK::wrap(p_kerneli->gram());
      break;
    default:
      return false;
      break;
  }
  if(p_kerneld) delete p_kerneld;
  if(p_kerneli) delete p_kerneli;
  // return result
  return true;
}


