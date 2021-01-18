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

/** @file STK_HDGaussian_ajk_bk_qk_d.h
 *  @brief In this file we define the HDGaussian_AjkBkQkDk class
 **/

#ifndef STK_HDGAUSSIAN_AJKBKQKD_H
#define STK_HDGAUSSIAN_AJKBKQKD_H

#include "../HDGaussianModels/STK_HDGaussianBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class HDGaussian_AjkBkQkD;

namespace hidden
{
/** @ingroup Clustering
 *  Traits class for the MixtureHDGaussian_ajk_bk_qk_d traits policy. */
template<class Array_>
struct MixtureTraits< HDGaussian_AjkBkQkD_<Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the parameters of a MixturGaussian_sjk model*/
  typedef ModelParameters<Clust::HDGaussian_AjkBkQkD_> Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a MixturGaussian_sj model.
 */
template<>
struct ModelParameters<Clust::HDGaussian_AjkBkQkD_>
{
    /** array of size nbCluster with the parameters mean of the variables */
    Array1D<CPointX> mean_;
    /** standard deviation of the variables */
    Array1D<CPointX>  sigma_;
    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster): mean_(nbCluster), sigma_(nbCluster) {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param)
                   : mean_(param.mean_), sigma_(param.sigma_)
    {}
    /** destructor */
    ~ModelParameters() {}
    /** @return the mean of the kth cluster and jth variable */
    inline Real const& mean(int k, int j) const { return mean_[k][j];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& sigma(int k, int j) const { return sigma_[k][j];}
    /** resize the set of parameter */
    inline void resize(Range const& range)
    {
      for (int k = mean_.begin(); k< mean_.end(); ++k)
      { mean_[k].resize(range) =0;
        sigma_[k].resize(range) = 1.;
      }
    }
};

/** @ingroup Clustering
 *  The diagonal Gaussian mixture model @c HDGaussian_AjkBkQkD is
 *  the most general HD Gaussian model and have a density function of the
 *  form
 **/
template<class Array>
class HDGaussian_AjkBkQkD: public HDGaussianBase<HDGaussian_AjkBkQkD<Array> >
{
  public:
    typedef DiagGaussianBase<HDGaussian_AjkBkQkDk<Array> > Base;
    using Base::param_;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    HDGaussian_AjkBkQkD( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    HDGaussian_AjkBkQkD( HDGaussian_AjkBkQkD const& model): Base(model) {}
    /** destructor */
    ~HDGaussian_AjkBkQkD() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0.;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      {
        if (param_.sigma_[k][j])
        { sum += Law::Normal::lpdf(p_data()->elt(i,j), param_.mean_[k][j], param_.sigma_[k][j]);}
      }
      return sum;
    }
    /** Initialize randomly the parameters of the Gaussian mixture. The centers
     *  will be selected randomly among the data set and the standard-deviation
     *  will be set to 1.
     */
    void randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** Compute the weighted mean and the common standard deviation. */
    bool run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return 2*this->nbCluster()*p_data()->sizeCols();}
};

/* Initialize randomly the parameters of the Gaussian mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void HDGaussian_AjkBkQkD<Array>::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{
  this->randomMean(p_tik);
  // compute the standard deviation
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    param_.sigma_[k] = Stat::varianceWithFixedMean(*p_data(), p_tik->col(k), param_.mean_[k], false).sqrt();
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("MixtureHDGaussian_ajk_bk_qk_d<Array>::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk)  done\n");
#endif
}

/* Compute the weighted means and the weighted standard deviations. */
template<class Array>
bool HDGaussian_AjkBkQkD<Array>::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{
  // compute the means
  if (!this->updateMean(p_tik)) return false;
  // compute the standard deviation
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    param_.sigma_[k] = Stat::varianceWithFixedMean(*p_data(), p_tik->col(k), param_.mean_[k], false).sqrt();
#ifdef STK_MIXTURE_DEBUG
    if( (param_.sigma_[k] <= 0).any()  )
    {
      stk_cout << _T("MixtureHDGaussian_ajk_bk_qk_d::run() failed\n");
      stk_cout << _T("p_tik->col(") << k << _T(") =\n") << p_tik->col(k).transpose() << _T("\n");
      stk_cout << _T("param_.mean_[") << k << _T("]  =") << param_.mean_[k];
      stk_cout << _T("param_.sigma_[") << k << _T("] =") << param_.sigma_[k];
    }
#endif
  }
  return true;
}

} // namespace STK

#endif /* STK_HDGAUSSIAN_AJKBKQKD_H */
