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

/** @file STK_DiagGaussian_sk.h
 *  @brief In this file we define the DiagGaussian_sk model
 **/

#ifndef STK_DIAGGAUSSIAN_SK_H
#define STK_DIAGGAUSSIAN_SK_H

#include "../DiagGaussianModels/STK_DiagGaussianBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class DiagGaussian_sk;

namespace hidden
{
/** @ingroup Clustering
 *  Traits class for the DiagGaussian_sk traits policy. */
template<class Array_>
struct MixtureTraits< DiagGaussian_sk<Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the parameters of a MixturGaussian_sk model*/
  typedef ModelParameters<Clust::Gaussian_sk_> Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  The diagonal Gaussian mixture model @c DiagGaussian_sk assumes an equal standard
 *  deviation in each cluster and has a density function of the form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d
 *    \frac{1}{\sqrt{2\pi}\sigma_{k}} \exp\left\{-\frac{(x^j-\mu^j_k)^2}{2(\sigma_{k})^2}\right\}.
 * \f]
 **/
template<class Array>
class DiagGaussian_sk: public DiagGaussianBase<DiagGaussian_sk<Array> >
{
  public:
    typedef DiagGaussianBase<DiagGaussian_sk<Array> > Base;
    using Base::param_;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    DiagGaussian_sk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    DiagGaussian_sk( DiagGaussian_sk const& model): Base(model) {}
    /** destructor */
    ~DiagGaussian_sk() {}
    /** Initialize randomly the parameters of the Gaussian mixture. The centers
     *  will be selected randomly among the data set and the standard-deviations
     *  will be set to 1.
     */
    void randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** Compute the weighted mean and the common standard deviation. */
    bool run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*p_data()->sizeCols() + this->nbCluster();}
};

/* Initialize randomly the parameters of the Gaussian mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void DiagGaussian_sk<Array>::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) 
{
  this->randomMean(p_tik);
  // compute the standard deviation
  Real variance;
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    variance = sqrt( ( p_tik->col(k).transpose()
                     *(*p_data() - (Const::Vector<Real>(this->nbSample()) * param_.mean_[k])
                      ).square()
                     ).sum() / (p_data()->sizeCols()*p_tk->elt(k))
                   );
    param_.sigma_[k] = ((variance<=0) || !Arithmetic<Real>::isFinite(variance))
                       ? 1.
                       : std::sqrt(variance/(this->nbSample()*p_data()->sizeCols()));
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("DiagGaussian_sk<Array>::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk)  done\n");
#endif
}

/* Compute the weighted mean and the common standard deviation. */
template<class Array>
bool DiagGaussian_sk<Array>::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) 
{
  // compute the means
  if (!this->updateMean(p_tik)) return false;
  // compute the standard deviation
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    param_.sigma_[k]
    = sqrt( ( p_tik->col(k).transpose()
             *(*p_data() - (Const::Vector<Real>(this->nbSample()) * param_.mean_[k])
              ).square()
            ).sum()
           /(p_data()->sizeCols()*p_tk->elt(k))
          );
//    if (param(k).sigma_ <= 0.) return false;
  }
  return true;
}

} // namespace STK

#endif /* STK_DiagGaussian_SK_H */
