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
 * created on: Oct 24, 2014
 * Author:   Serge Iovleff
 **/

/** @file STK_KernelParameters.cpp
 *  @brief In this file we implement the ModelParameters class for kernel
 *  mixture models
 **/


#include <Clustering/include/KernelModels/STK_KernelParameters.h>

namespace STK
{
/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Kmm_s_>::ModelParameters(int nbCluster)
               : sigma2_(1.), dim_(nbCluster, 10)
               , stat_sigma2_(), stat_dim_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Kmm_s_>::ModelParameters( ModelParameters const& param)
               : sigma2_(param.sigma2_), dim_(param.dim_)
               , stat_sigma2_(param.stat_sigma2_)
               , stat_dim_(param.stat_dim_)
{}
/* destructor */
ModelParameters<Clust::Kmm_s_>::~ModelParameters() {}

/* update statistics of the parameters. */
void ModelParameters<Clust::Kmm_s_>::updateStatistics()
{
  stat_sigma2_.update(sigma2_);
  for(int k=stat_dim_.begin(); k<stat_dim_.end(); ++k)
  { stat_dim_[k].update(dim_[k]);}
}
/* Set the computed statistics */
void ModelParameters<Clust::Kmm_s_>::setStatistics()
{
  sigma2_ = stat_sigma2_.mean();
  stat_sigma2_.release();
  for(int k=stat_dim_.begin(); k<stat_dim_.end(); ++k)
  {
    dim_[k] = stat_dim_[k].mean();
    stat_dim_[k].release();
  }
}
/* Release the computed statistics */
void ModelParameters<Clust::Kmm_s_>::releaseStatistics()
{
  stat_sigma2_.release();
  for(int k=stat_dim_.begin(); k<stat_dim_.end(); ++k)
  { stat_dim_[k].release();}
}


/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Kmm_sk_>::ModelParameters(int nbCluster)
               : sigma2_(nbCluster, 1.), dim_(nbCluster, 10)
               , stat_sigma2_(nbCluster), stat_dim_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Kmm_sk_>::ModelParameters( ModelParameters const& param)
               : sigma2_(param.sigma2_), dim_(param.dim_)
               , stat_sigma2_(param.stat_sigma2_)
               , stat_dim_(param.stat_dim_)
{}
/* destructor */
ModelParameters<Clust::Kmm_sk_>::~ModelParameters() {}

/* update statistics of the parameters. */
void ModelParameters<Clust::Kmm_sk_>::updateStatistics()
{
  for(int k=stat_dim_.begin(); k<stat_dim_.end(); ++k)
  {
    stat_sigma2_[k].update(sigma2_[k]);
    stat_dim_[k].update(dim_[k]);}
}
/* Set the computed statistics */
void ModelParameters<Clust::Kmm_sk_>::setStatistics()
{
  for(int k=stat_dim_.begin(); k<stat_dim_.end(); ++k)
  {
    sigma2_[k] = stat_sigma2_[k].mean();
    stat_sigma2_[k].release();
    dim_[k] = stat_dim_[k].mean();
    stat_dim_[k].release();
  }
}
/* Release the computed statistics */
void ModelParameters<Clust::Kmm_sk_>::releaseStatistics()
{
  for(int k=stat_dim_.begin(); k<stat_dim_.end(); ++k)
  {
    stat_sigma2_[k].release();
    stat_dim_[k].release();
  }
}


} // namespace STK
