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

/** @file STK_Categorical_pjk.h
 *  @brief In this file we implement the Categorical_pjk mixture model
 **/

#ifndef STK_CATEGORICAL_PJK_H
#define STK_CATEGORICAL_PJK_H

#include "../CategoricalModels/STK_CategoricalBase.h"

namespace STK
{

//forward declaration
template<class Array>class Categorical_pjk;

namespace hidden
{
/** @ingroup Clustering
 *  Traits class for the Categorical_pjk traits policy. */
template<class Array_>
struct MixtureTraits< Categorical_pjk<Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the parameters of a Categorical_pjk model*/
  typedef ModelParameters<Clust::Categorical_pjk_> Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  The diagonal Categorical mixture model @c Categorical_pjk is
 *  the most general diagonal Categorical model and has a probability
 *  function of the form
 * \f[
 *    P(\mathbf{x}=(l_1,\ldots,l_d)|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d p_{kl_j}^j.
 * \f]
 **/
template<class Array>
class Categorical_pjk: public CategoricalBase<Categorical_pjk<Array> >
{
  public:
    typedef CategoricalBase<Categorical_pjk<Array> > Base;
    using Base::param_;
    using Base::p_data;
    using Base::modalities_;

    /** default constructor
     * @param nbCluster number of cluster of the model
     **/
    Categorical_pjk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model model to copy
     **/
    Categorical_pjk( Categorical_pjk const& model): Base(model) {}
    /** destructor */
    ~Categorical_pjk() {}
//    Real lnComponentProbability(int i, int k) const;
    /** Initialize randomly the parameters of the Categorical mixture. */
    void randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** Compute the weighted probabilities. */
    bool run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*((this->nbModalities_-1).sum());}
};

/* Initialize randomly the parameters of the Categorical mixture. */
template<class Array>
void Categorical_pjk<Array>::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{
  for (int k = p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    param_.proba_[k].randUnif();
    for (int j=param_.proba_[k].beginCols(); j< param_.proba_[k].endCols(); ++j)
    { param_.proba_[k].col(j) /= param_.proba_[k].col(j).sum();}
  }
}


/* Compute the modalities probabilities */
template<class Array>
bool Categorical_pjk<Array>::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{
  for (int k = p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    param_.proba_[k] = 0.;
    for (int j = p_data()->beginCols(); j < p_data()->endCols(); ++j)
    {
      // count the number of modalities weighted by the tik
      for (int i = p_data()->beginRows(); i < p_data()->endRows(); ++i)
      { param_.proba_[k](p_data()->elt(i, j), j) += p_tik->elt(i, k);}
      // normalize the probabilities
      Real sum = param_.proba_[k].col(j).sum();
      if (sum) { param_.proba_[k].col(j) /= sum;}
    }
  }
  return true;
}

} // namespace STK

#endif /* STK_Categorical_PJK_H */
