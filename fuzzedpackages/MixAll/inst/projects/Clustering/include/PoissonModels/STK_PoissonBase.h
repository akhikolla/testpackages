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

/* Project:  stkpp::Clustering
 * created on: Dec 9, 2014
 * Authors: Serge Iovleff
 **/

/** @file STK_PoissonBase.h
 *  @brief In this file we implement the base class for the poisson mixture models
 **/

#ifndef STK_POISSONBASE_H
#define STK_POISSONBASE_H

#include "../STK_IMixtureDensity.h"

#include <STatistiK/include/STK_Law_Poisson.h>
#include "../PoissonModels/STK_PoissonParameters.h"

namespace STK
{
/** @ingroup Clustering
 *  Base class for the Poisson models
 **/
template<class Derived>
class PoissonBase: public IMixtureDensity<Derived >
{
  public:
    typedef IMixtureDensity<Derived > Base;
    using Base::param_;
    using Base::p_data;
    using Base::nbCluster;

  protected:
    /** default constructor
     *  @param nbCluster number of cluster in the model
     **/
    PoissonBase( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    PoissonBase( PoissonBase const& model): Base(model) {}
    /** destructor */
    ~PoissonBase() {}

  public:
    /** @return the value of lambda of the kth cluster and jth variable */
    inline Real lambda(int k, int j) const { return param_.lambda(k,j);}
    /** Initialize the parameters of the model. */
    void initializeModelImpl() { param_.resize(p_data()->cols());}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    Real lnComponentProbability(int i, int k) const;
    /** @return an imputation value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute
     *  @param pk the probabilities of each class for the ith individual
     **/
    template<class Weights>
    int impute(int i, int j, Weights const& pk) const;
    /** @return a value to impute for the jth variable of the ith sample*/
    Real impute(int i, int j, CArrayXX const*  p_tik) const;
    /** @return a simulated value for the jth variable of the ith sample
     *  in the kth cluster.
     *  @param i,j,k indexes of the data to simulate
     **/
    inline int rand(int i, int j, int k) const
    { return Law::Poisson::rand(lambda(k,j));}
    /** This function is used in order to get the current values of the lambdas.
     *  @param[out] params the array with the parameters of the mixture.
     */
    template<class Array>
    void getParameters(Array& params) const;
    /** This function can be used to write summary of parameters to the output stream.
     *  @param p_tik a constant pointer on the posterior probabilities
     *  @param os Stream where you want to write the summary of parameters.
     */
    void writeParameters(CArrayXX const* p_tik, ostream& os) const;
};

/* @return the value of the probability of the i-th sample in the k-th component.
 *  @param i,k indexes of the sample and of the component
 **/
template<class Derived>
Real PoissonBase<Derived>::lnComponentProbability(int i, int k) const
{
  Real sum =0.;
  for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
  {
    Real value = lambda(k,j);
    if (value)
    { sum += Law::Poisson::lpdf(p_data()->elt(i,j), value);}
  }
  return sum;
}
/* Implementation  */
template<class Derived>
template<class Weights>
int PoissonBase<Derived>::impute(int i, int j, Weights const& pk) const
{
  Real sum = 0.;
  for (int k= pk.begin(); k < pk.end(); ++k)
  { sum += pk[k] * lambda(k,j);}
  return std::floor(sum+0.5);
}

/* @return a value to impute for the jth variable of the ith sample*/
template<class Derived>
Real PoissonBase<Derived>::impute(int i, int j, CArrayXX const*  p_tik) const
{
  Real sum = 0.;
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  { sum += p_tik->elt(i,k) * lambda(k,j);}
  return sum;
}

/* This function is used in order to get the current values of the lambdas.
 *  @param[out] params the array with the parameters of the mixture.
 */
template<class Derived>
template<class Array>
void PoissonBase<Derived>::getParameters(Array& params) const
{
  params.resize(nbCluster(), p_data()->cols());
  for (int k= params.beginRows(); k < params.endRows(); ++k)
  {
    for (int j= p_data()->beginCols();  j < p_data()->endCols(); ++j)
    { params(k, j) = lambda(k,j);}
  }
}


/* This function can be used to write summary of parameters to the output stream.
 *  @param p_tik a constant pointer on the posterior probabilities
 *  @param os Stream where you want to write the summary of parameters.
 */
template<class Derived>
void PoissonBase<Derived>::writeParameters(CArrayXX const* p_tik, ostream& os) const
{
  CPointX params(p_data()->cols());
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    for (int j= p_data()->beginCols();  j < p_data()->endCols(); ++j)
    { params[j] = lambda(k,j);}
    os << _T("---> Component ") << k << _T("\n");
    os << _T("lambda = ") << params;
  }
}


} // namespace STK

#endif /* STK_PoissonBASE_H */
