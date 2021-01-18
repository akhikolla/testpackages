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

/** @file STK_Poisson_ljk.h
 *  @brief In this file we implement the Poisson_ljk class
 **/

#ifndef STK_POISSON_LJK_H
#define STK_POISSON_LJK_H

#include "../PoissonModels/STK_PoissonBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class Poisson_ljk;

namespace hidden
{
/** @ingroup Clustering
 *  Traits class for the Poisson_ljk traits policy. */
template<class Array_>
struct MixtureTraits< Poisson_ljk<Array_> >
{
  typedef Array_                  Array;
  typedef typename Array::Type    Type;
  /** Type of the structure storing the parameters of a MixturePoisson_ljk model*/
  typedef ModelParameters<Clust::Poisson_ljk_> Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  The Poisson mixture model @c Poisson_ljk is the most general Poisson model
 *  and have a probability function of the form
 * \f[
 *    P(\mathbf{x}=(n_1,\ldots,n_d)|\theta)
 *     = \sum_{k=1}^K p_k \prod_{j=1}^d e^{-\lambda_{jk}} \frac{\lambda_{jk}^{n_j}}{n_j!}.
 * \f]
 **/
template<class Array>
class Poisson_ljk: public PoissonBase<Poisson_ljk<Array> >
{
  public:
    typedef PoissonBase<Poisson_ljk<Array> > Base;
    typedef typename hidden::MixtureTraits< Poisson_ljk<Array> >::Parameters Parameters;
    using Base::param_;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Poisson_ljk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    Poisson_ljk( Poisson_ljk const& model): Base(model) {}
    /** destructor */
    ~Poisson_ljk() {}
    /** Initialize randomly the parameters of the Poisson mixture. */
    void randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** Compute the MStep for the lambda_jk
     *  \f$ \lambda_jk = \sum_i X_{ij} t_{ik} \f$ */
    bool run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*p_data()->sizeCols();}
};

/* Initialize randomly the parameters of the Poisson mixture. */
template<class Array>
void Poisson_ljk<Array>::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{
  for (int j=p_data()->beginCols(); j< p_data()->endCols(); ++j)
  {
    Real m = p_data()->col(j).template cast<Real>().mean();
    for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
    { param_.lambda_[k][j] = Law::Exponential::rand(m);}
  }
}


/* Compute the lambda_jk  */
template<class Array>
bool Poisson_ljk<Array>::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{
  for (int j=p_data()->beginCols(); j< p_data()->endCols(); ++j)
  {
    for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
    {
      param_.lambda_[k][j] = 0.;
      for (int i= p_tik->beginRows(); i < p_tik->endRows(); ++i)
      { param_.lambda_[k][j] += p_data()->elt(i,j) * p_tik->elt(i,k);}
      param_.lambda_[k][j] /= p_tk->elt(k);
    }
  }
  return true;
}

} // namespace STK

#endif /* STK_Poisson_LJK_H */
