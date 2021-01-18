/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016 Serge Iovleff

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: Oct 24, 2014
 * Author:   Serge Iovleff
 **/

/** @file STK_Kmm_s.cpp
 *  @brief In this file we implement the Kmm_s class
 **/


#include <Analysis/include/STK_Const_Math.h>
#include <Clustering/include/KernelModels/STK_Kmm_s.h>
#include <STatistiK/include/STK_Law_Normal.h>
#include <STatistiK/include/STK_Stat_Functors.h>

namespace STK
{
/* default constructor
 * @param nbCluster number of cluster in the model
 **/
Kmm_s::Kmm_s( int nbCluster): Base(nbCluster) {}
/* copy constructor
 *  @param model The model to copy
 **/
Kmm_s::Kmm_s( Kmm_s const& model): Base(model) {}
/* destructor */
Kmm_s::~Kmm_s() {}
/* @return the number of free parameters of the model */
int Kmm_s::computeNbFreeParameters() const
{ return param_.dim_.sum() + 1;}

/* @return the value of the probability of the i-th sample in the k-th component.
 *  @param i,k indexes of the sample and of the component
 **/
Real Kmm_s::lnComponentProbability(int i, int k) const
{
  return(- dik_.elt(i,k)/(2.*param_.sigma2_)
         - (std::log(param_.sigma2_)+2.*Const::_LNSQRT2PI_)*param_.dim_[k]/2.);
}

/* Initialize randomly the parameters of the Gaussian mixture. */
void Kmm_s::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{
#if STK_Kernel_DEBUG | STK_MIXTURE_VERBOSE
  stk_cout << _T("Entering Kmm_s::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk)\n");
#endif
  // compute the standard deviation
  compute_dik(p_tik, p_tk);
  param_.sigma2_ = dik_.prod(*p_tik).sum()/(this->nbSample() * param_.dim_.sum())
                 + std::abs(Law::generator.randGauss(0, 0.05));
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("Kmm_s::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk)  done\n");
  stk_cout << param_.sigma2_ << "\n";
#endif
}

/* Compute the weighted means and the weighted standard deviations. */
bool Kmm_s::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{
#if STK_Kernel_DEBUG | STK_MIXTURE_VERBOSE
  stk_cout << _T("Entering Kmm_s::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk)\n");
#endif
  compute_dik(p_tik, p_tk);
  param_.sigma2_ =  p_tik->prod(dik_).sum()/p_tk->dot(param_.dim_);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("Kmm_s::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) done\n");
  stk_cout << param_.sigma2_ << "\n";
#endif
  return (param_.sigma2_ <= 0.) ? false : true;
}

} // namespace STK
