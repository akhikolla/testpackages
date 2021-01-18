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

/** @file STK_Gamma_ajk_bj.h
 *  @brief In this file we define the Gamma_pk_ajk_bj and Gamma_p_ajk_bj mixture models.
 **/

#ifndef STK_GAMMA_AJK_BJ_H
#define STK_GAMMA_AJK_BJ_H

#include <STatistiK/include/STK_Law_Exponential.h>
#include "../GammaModels/STK_GammaBase.h"

#define MAXITER 400
#define TOL 1e-8

namespace STK
{
template<class Array>class Gamma_ajk_bj;

namespace hidden
{
/** @ingroup Clustering
 *  Traits class for the Gamma_ajk_bj traits policy. */
template<class Array_>
struct MixtureTraits< Gamma_ajk_bj<Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the parameters of a Gamma_ajk_bj model*/
  typedef ModelParameters<Clust::Gamma_ajk_bj_> Parameters;
};

} // namespace Clust

/** @ingroup Clustering
 *  Gamma_ajk_bj is a mixture model of the following form
 * \f[
 *     f(\mathbf{x}_i|\theta) = \sum_{k=1}^K p_k
 *     \prod_{j=1}^p\left(\frac{x_i^j}{b_{j}}\right)^{a_{jk}-1}
 *                   \frac{e^{-x_i^j/b_{j}}} {b_{j} \, \Gamma(a_{jk})},
 *      \quad x_i^j>0, \quad i=1,\ldots,n.
 * \f]
 **/
template<class Array>
class Gamma_ajk_bj: public GammaBase<Gamma_ajk_bj<Array> >
{
  public:
    typedef GammaBase<Gamma_ajk_bj<Array> > Base;
    using Base::param_;
    
    using Base::p_data;
    using Base::meanjk;
    using Base::variancejk;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Gamma_ajk_bj( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    Gamma_ajk_bj( Gamma_ajk_bj const& model): Base(model) {}
    /** destructor */
    ~Gamma_ajk_bj() {}
    /** Initialize randomly the parameters of the Gamma mixture. */
    void randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** Compute the weighted mean and the common variance. */
    bool run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*p_data()->sizeCols()+ p_data()->sizeCols();}
};

/* Initialize randomly the parameters of the gamma mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void Gamma_ajk_bj<Array>::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) 
{
    // compute moments
    this->moments(p_tik);
  for (int j=p_data()->beginCols(); j < p_data()->endCols(); ++j)
  {
    Real value =0.;
    for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
    {
      Real mean = meanjk(j,k), variance = variancejk(j,k);
      param_.shape_[k][j] = Law::Exponential::rand((mean*mean/variance));
      value += p_tk->elt(k) * variance/mean;
    }
    param_.scale_[j] = Law::Exponential::rand(value/(this->nbSample()));
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T(" Gamma_ajk_bj<Array>::randomInit done\n");
#endif
}

/* Compute the weighted mean and the common variance. */
template<class Array>
bool Gamma_ajk_bj<Array>::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) 
{
  if (!this->moments(p_tik)) { return false;}
  // start estimations of the ajk and bj
  Real qvalue = this->qValue(p_tik, p_tk);
  int iter;
  for(iter=0; iter<MAXITER; ++iter)
  {
    for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
    {
      Real num=0., den = 0.;
      // compute ajk
      for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
      {
        // moment estimate and oldest value
        Real x0 = this->meanjk(j,k)*this->meanjk(j,k)/this->variancejk(j,k);
        Real x1 = param_.shape_[k][j];
        if ((x0 <=0.) || !Arithmetic<Real>::isFinite(x0)) return false;
        // compute shape
        hidden::invPsi f(param_.meanLog_[k][j] - std::log(param_.scale_[j]));
        Real a =  Algo::findZero(f, x0, x1, TOL);

        if (!Arithmetic<Real>::isFinite(a))
        {
          param_.shape_[k][j] = x0; // use moment estimate
#ifdef STK_MIXTURE_DEBUG
          stk_cout << _T("ML estimation failed in Gamma_ajk_bj::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) \n");
          stk_cout << "x0 =" << x0 << _T("\n";);
          stk_cout << "f(x0) =" << f(x0) << _T("\n";);
          stk_cout << "x1 =" << x1 << _T("\n";);
          stk_cout << "f(x1) =" << f(x1) << _T("\n";);
#endif
        }
        else { param_.shape_[k][j] = a;}
        num += param_.mean_[k][j]  * p_tk->elt(k);
        den += param_.shape_[k][j] * p_tk->elt(k);
      }
      // compute b_j
      Real b = num/den;
      // divergence
      if (!Arithmetic<Real>::isFinite(b)) { return false;}
      param_.scale_[j] = b;
    }
    // check convergence
    Real value = this->qValue(p_tik, p_tk);
#ifdef STK_MIXTURE_DEBUG
    if (value < qvalue)
    {
      stk_cout << _T("In Gamma_ajk_bj::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) : run( CArrayXX const* const& p_tik, CPointX const* const& p_tk)  diverge\n");
      stk_cout << _T("New value =") << value << _T(", qvalue =") << qvalue << _T("\n");
    }
#endif
    if ((value - qvalue) < TOL) break;
    qvalue = value;
  }
#ifdef STK_MIXTURE_DEBUG
  if (iter == MAXITER)
  {
    stk_cout << _T("In Gamma_ajk_bj::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) : run( CArrayXX const* const& p_tik, CPointX const* const& p_tk)  did not converge\n");
    stk_cout << _T("qvalue =") << qvalue << _T("\n");
  }
#endif
  return true;
}

}  // namespace STK

#undef MAXITER
#undef TOL

#endif /* STK_Gamma_AJK_BJ_H */
