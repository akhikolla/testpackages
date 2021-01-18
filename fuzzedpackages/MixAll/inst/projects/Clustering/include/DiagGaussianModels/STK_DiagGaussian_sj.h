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

/** @file STK_DiagGaussian_sj.h
 *  @brief In this file we define and implement the DiagGaussian_sj class
 **/

#ifndef STK_DIAGGAUSSIAN_SJ_H
#define STK_DIAGGAUSSIAN_SJ_H

#include "../DiagGaussianModels/STK_DiagGaussianBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class DiagGaussian_sj;

namespace hidden
{
/** @ingroup hidden
 *  Traits class for the Gaussian_s traits policy. */
template<class Array_>
struct MixtureTraits< DiagGaussian_sj<Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the parameters of a MixturGaussian_sj model*/
  typedef ModelParameters<Clust::Gaussian_sj_> Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  The diagonal Gaussian mixture model DiagGaussian_sj have a density function of the
 *  form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d
 *    \frac{1}{\sqrt{2\pi}\sigma_j} \exp\left\{-\frac{(x^j-\mu^j_k)^2}{2\sigma_j^2}\right\}.
 * \f]
 **/
template<class Array>
class DiagGaussian_sj: public DiagGaussianBase<DiagGaussian_sj<Array> >
{
  public:
    typedef DiagGaussianBase<DiagGaussian_sj<Array> > Base;
    using Base::param_;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    DiagGaussian_sj( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    DiagGaussian_sj( DiagGaussian_sj const& model)
                     : Base(model)
    {}
    /** destructor */
    ~DiagGaussian_sj() {}
    /** Initialize randomly the parameters of the Gaussian mixture. The centers
     *  will be selected randomly among the data set and the standard-deviation
     *  will be set to 1.
     */
    void randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** Compute the weighted mean and the common standard deviation. */
    bool run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*p_data()->sizeCols()+p_data()->sizeCols();}
};

/* Initialize randomly the parameters of the Gaussian mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void DiagGaussian_sj<Array>::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) 
{
  // compute the initial mean
  this->randomMean(p_tik);
  // compute the standard deviation
  Array2DPoint<Real> variance(p_data()->cols(), 0.);
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    variance += p_tik->col(k).transpose()
               *(*p_data() - (Const::Vector<Real>(this->nbSample()) * param_.mean_[k])
                ).square()
                ;
  }
  // store the standard deviation
  param_.sigma_ = (variance /= this->nbSample()).sqrt();
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("DiagGaussian_sj<Array>::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk)  done\n");
#endif
}

/* Compute the weighted mean and the common standard deviation. */
template<class Array>
bool DiagGaussian_sj<Array>::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) 
{
  // compute the means
  if (!this->updateMean(p_tik)) return false;
  // compute the standard deviation
  Array2DPoint<Real> variance(p_data()->cols(), 0.);
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    variance += p_tik->col(k).transpose()
               *(*p_data() - (Const::Vector<Real>(this->nbSample()) * param_.mean_[k])
                ).square()
                ;
  }
//  if (variance.nbAvailableValues() != p_data()->sizeCols()) return false;
//  if ((variance > 0.).template cast<int>().sum() != p_data()->sizeCols()) return false;
  // compute the standard deviation
  param_.sigma_ = (variance /= this->nbSample()).sqrt();
  return true;
}

} // namespace STK

#endif /* STK_DiagGaussian_SJ_H */
