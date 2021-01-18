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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project: stkpp::Clustering
 * created on: 5 sept. 2013
 * Author:  iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Gamma_ak_bjk.h
 *  @brief In this file we define the Gamma_pk_ak_bjk and Gamma_p_ak_bjk models.
 **/

#ifndef STK_GAMMA_AK_BJK_H
#define STK_GAMMA_AK_BJK_H

#include <STatistiK/include/STK_Law_Exponential.h>
#include "../GammaModels/STK_GammaBase.h"

namespace STK
{
template<class Array>class Gamma_ak_bjk;

namespace hidden
{
/** @ingroup Clustering
 * Traits class for the Gamma_ak_bjk traits policy
 **/
template<class Array_>
struct MixtureTraits< Gamma_ak_bjk<Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the parameters of a Gamma_ak_bjk model*/
  typedef ModelParameters<Clust::Gamma_ak_bjk_> Parameters;
};

} // namespace Clust

/** @ingroup Clustering
 *  Gamma_ak_bjk is a mixture model of the following form
 * \f[
 *     f(\mathbf{x}_i|\theta) = \sum_{k=1}^K p_k
 *     \prod_{j=1}^p\left(\frac{x_i^j}{b_{jk}}\right)^{a_{k}-1}
 *                   \frac{e^{-x_i^j/b_{jk}}}{b_{jk} \, \Gamma(a_{k})},
 *      \quad x_i^j>0, \quad i=1,\ldots,n.
 * \f]
 **/
template<class Array>
class Gamma_ak_bjk: public GammaBase< Gamma_ak_bjk<Array> >
{
  public:
    typedef GammaBase< Gamma_ak_bjk<Array> > Base;
    using Base::param_;
    using Base::p_data;
    using Base::meanjk;
    using Base::variancejk;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Gamma_ak_bjk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    Gamma_ak_bjk( Gamma_ak_bjk const& model): Base(model) {}
    /** destructor */
    ~Gamma_ak_bjk() {}
    /** Initialize randomly the parameters of the Gamma mixture. The shape
     *  will be selected randomly using an exponential of parameter mean^2/variance
     *  and the scale will be selected randomly using an exponential of parameter
     *  variance/mean.
     */
    void randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** Compute the run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) . */
    bool run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*p_data()->sizeCols()+this->nbCluster();}
};

/* Initialize randomly the parameters of the Gaussian mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void Gamma_ak_bjk<Array>::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) 
{
    // compute moments
    this->moments(p_tik);
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    Real value = 0.0;
    for (int j=p_data()->beginCols(); j < p_data()->endCols(); ++j)
    {
      Real mean = meanjk(j,k), variance = variancejk(j,k);
      param_.scale_[k][j] = Law::Exponential::rand((variance/mean));
      value += mean*mean/variance;
    }
    param_.shape_[k]= Law::Exponential::rand(value/(p_data()->sizeCols()));
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T(" Gamma_ak_bjk<Array>::randomInit done\n");
#endif
}

/* Compute the weighted mean and the common variance. */
template<class Array>
bool Gamma_ak_bjk<Array>::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) 
{
  if (!this->moments(p_tik)) { return false;}
  // estimate a and b
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    // moment estimate and oldest value
    Real x0 = (param_.mean_[k].square()/param_.variance_[k]).mean();
    Real x1 =  param_.shape_[k];
    if ((x0 <=0.) || (isNA(x0))) return false;

    // get shape
    hidden::invPsiMLog f( (param_.meanLog_[k]-param_.mean_[k].log()).mean() );
    Real a = Algo::findZero(f, x0, x1, 1e-08);
    if (!Arithmetic<Real>::isFinite(a))
    {
#ifdef STK_MIXTURE_DEBUG
        stk_cout << "ML estimation failed in Gamma_ak_bjk::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) \n";
        stk_cout << "x0 =" << x0 << _T("\n";);
        stk_cout << "f(x0) =" << f(x0) << _T("\n";);
        stk_cout << "x1 =" << x1 << _T("\n";);
        stk_cout << "f(x1) =" << f(x1) << _T("\n";);
#endif
      a = x0; // use moment estimate
    }
    // set values
    param_.shape_[k]= a;
    for (int j=p_data()->beginCols(); j < p_data()->endCols(); ++j)
    { param_.scale_[k][j] = param_.mean_[k][j]/a; }
  }
  return true;
}

}  // namespace STK

#endif /* STK_Gamma_AK_BJK_H */
