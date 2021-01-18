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

/** @file  STK_PoissonParameters.h
 *  @brief In this file we implement the Parameters classes for Poisson
 *  mixture models
 **/

#include "../include/PoissonModels/STK_PoissonParameters.h"

namespace STK
{
/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Poisson_ljk_>::ModelParameters( int nbCluster)
                                                     : lambda_(nbCluster)
                                                     , stat_lambda_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Poisson_ljk_>::ModelParameters( ModelParameters const& param)
                                                     : lambda_(param.lambda_)
                                                     , stat_lambda_(param.stat_lambda_)
{}
/* destructor */
ModelParameters<Clust::Poisson_ljk_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Poisson_ljk_>::resize(Range const& range)
{
  for (int k = lambda_.begin(); k< lambda_.end(); ++k)
  {
    lambda_[k].resize(range) = 1.;
    stat_lambda_[k].resize(range);
  }
}
/* update statistics of the parameters. */
void ModelParameters<Clust::Poisson_ljk_>::updateStatistics()
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  { stat_lambda_[k].update(lambda_[k]);}
}
/* Set the computed statistics */
void ModelParameters<Clust::Poisson_ljk_>::setStatistics()
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  {
    lambda_[k] = stat_lambda_[k].mean();
    stat_lambda_[k].release();
  }
}
/* Release the computed statistics */
void ModelParameters<Clust::Poisson_ljk_>::releaseStatistics()
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  { stat_lambda_[k].release();}
}


/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Poisson_ljlk_>::ModelParameters(int nbCluster)
                : lambdak_(nbCluster), lambdaj_()
                , stat_lambdak_(nbCluster), stat_lambdaj_()
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Poisson_ljlk_>::ModelParameters( ModelParameters const& param)
               : lambdak_(param.lambdak_)
               , lambdaj_(param.lambdaj_)
               , stat_lambdak_(param.stat_lambdak_)
               , stat_lambdaj_(param.stat_lambdaj_)
{}
/* destructor */
ModelParameters<Clust::Poisson_ljlk_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Poisson_ljlk_>::resize(Range const& range)
{
  for (int k = lambdak_.begin(); k< lambdak_.end(); ++k)
  {
    lambdak_[k] = 1.;
    stat_lambdak_[k].release();
  }
  lambdaj_.resize(range) = 1;
  stat_lambdaj_.resize(range);
}

/* update statistics of the parameters. */
void ModelParameters<Clust::Poisson_ljlk_>::updateStatistics()
{
  for(int k=stat_lambdak_.begin(); k<stat_lambdak_.end(); ++k)
  { stat_lambdak_[k].update(lambdak_[k]);}
  stat_lambdaj_.update(lambdaj_);
}
/* Set the computed statistics */
void ModelParameters<Clust::Poisson_ljlk_>::setStatistics()
{
  for(int k=stat_lambdak_.begin(); k<stat_lambdak_.end(); ++k)
  {
    lambdak_[k] = stat_lambdak_[k].mean();
    stat_lambdak_[k].release();
  }
  lambdaj_ = stat_lambdaj_.mean();
  stat_lambdaj_.release();

}
/* Release the computed statistics */
void ModelParameters<Clust::Poisson_ljlk_>::releaseStatistics()
{
  for(int k=stat_lambdak_.begin(); k<stat_lambdak_.end(); ++k)
  { stat_lambdak_[k].release();}
  stat_lambdaj_.release();
}

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Poisson_lk_>::ModelParameters(int nbCluster)
               : lambda_(nbCluster), stat_lambda_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Poisson_lk_>::ModelParameters( ModelParameters const& param)
               : lambda_(param.lambda_), stat_lambda_(param.stat_lambda_)
{}
/* destructor */
ModelParameters<Clust::Poisson_lk_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Poisson_lk_>::resize(Range const& range)
{
  for (int k = lambda_.begin(); k< lambda_.end(); ++k)
  { lambda_[k] = 1.;
   stat_lambda_[k].release();
  }
}

/* update statistics of the parameters. */
void ModelParameters<Clust::Poisson_lk_>::updateStatistics()
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  { stat_lambda_[k].update(lambda_[k]);}
}
/* Set the computed statistics */
void ModelParameters<Clust::Poisson_lk_>::setStatistics()
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  {
    lambda_[k] = stat_lambda_[k].mean();
    stat_lambda_[k].release();
  }
}
/* Release the computed statistics */
void ModelParameters<Clust::Poisson_lk_>::releaseStatistics()
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  { stat_lambda_[k].release();}
}

} // namespace STK

