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
 * created on: Oct 24, 2013
 * Author:   Serge Iovleff
 **/

/** @file  STK_DiagGaussian_s.h
 *  @brief In this file we implement the DiagGaussian_s class
 **/

#ifndef STK_DIAGGAUSSIAN_S_H
#define STK_DIAGGAUSSIAN_S_H

#include "../DiagGaussianModels/STK_DiagGaussianBase.h"

namespace STK
{

// forward declaration
template<class Array>class DiagGaussian_s;

namespace hidden
{
/** @ingroup hidden
 *  Traits class for the DiagGaussian_s traits policy. */
template<class Array_>
struct MixtureTraits< DiagGaussian_s<Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gaussian_s_> Parameters;
};

} // namespace Clust

/** @ingroup Clustering
 *  The diagonal DiagGaussian_s mixture model hse a density function of the form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d
 *    \frac{1}{\sqrt{2\pi}\sigma} \exp\left\{-\frac{(x^j-\mu^j_k)^2}{2\sigma^2}\right\}.
 * \f]
 **/
template<class Array>
class DiagGaussian_s: public DiagGaussianBase<DiagGaussian_s<Array> >
{
  public:
    typedef DiagGaussianBase<DiagGaussian_s<Array> > Base;
    using Base::param_;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    DiagGaussian_s( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    DiagGaussian_s( DiagGaussian_s const& model): Base(model) {}
    /** destructor */
    ~DiagGaussian_s() {}
    /** Initialize randomly the parameters of the Gaussian mixture. The centers
     *  will be selected randomly among the data set and the standard-deviation
     *  will be set to 1.
     */
    void randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** Compute the weighted mean and the common standard deviation. */
    bool run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*p_data()->sizeCols()+1;}
};

/* Initialize randomly the parameters of the Gaussian mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void DiagGaussian_s<Array>::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{
  this->randomMean(p_tik);
  // compute the standard deviation
  Real variance = 0.0;
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    variance += (  p_tik->col(k).transpose()
                 * (*p_data() - (Const::Vector<Real>(this->nbSample()) * param_.mean_[k])
                   ).square()
                ).sum();
  }
  param_.sigma_ = std::sqrt(variance/(this->nbSample()*p_data()->sizeCols()));
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Gaussian_s<Array>::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk)  done\n");
#endif
}

/* Compute the weighted mean and the common standard deviation. */
template<class Array>
bool DiagGaussian_s<Array>::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{
  // compute the means
  if (!this->updateMean(p_tik)) return false;
  // compute the standard deviation
  Real variance = 0.0;
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    variance += ( p_tik->col(k).transpose()
                 * (*p_data() - (Const::Vector<Real>(this->nbSample()) * param_.mean_[k])
                   ).square()
                ).sum();
  }
  param_.sigma_ = std::sqrt(variance/(this->nbSample()*p_data()->sizeCols()));
#ifdef STK_MIXTURE_DEBUG
    if( param_.sigma_ <= 0  )
    {
      stk_cout << _T("DiagGaussian_s::run() failed\n");
      stk_cout << _T("param_.mean_ =") << param_.mean_;
      stk_cout << _T("param_.sigma_=") << param_.sigma_ << _T("\n");
    }
#endif
  //if ((variance<=0) || !Arithmetic<Real>::isFinite(variance)) return false;
  return true;
}

} // namespace STK

#endif /* STK_DiagGaussian_S_H */
