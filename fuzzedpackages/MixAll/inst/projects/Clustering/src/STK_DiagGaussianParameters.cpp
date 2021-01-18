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

/** @file  STK_DiagGaussianParameters.cpp
 *  @brief In this file we implement the Parameters classes for Diagonal Gaussian
 *  mixture models
 **/

#include "../include/DiagGaussianModels/STK_DiagGaussianParameters.h"

namespace STK
{
/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gaussian_sjk_>::ModelParameters(int nbCluster)
                                                      : mean_(nbCluster), sigma_(nbCluster)
                                                      , stat_mean_(nbCluster), stat_sigma_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gaussian_sjk_>::ModelParameters( ModelParameters const& param)
               : mean_(param.mean_), sigma_(param.sigma_)
               , stat_mean_(param.stat_mean_), stat_sigma_(param.stat_sigma_)
{}
/* destructor */
ModelParameters<Clust::Gaussian_sjk_>::~ModelParameters() {}
/* copy operator.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gaussian_sjk_>& ModelParameters<Clust::Gaussian_sjk_>::operator=( ModelParameters const& param)
{
  mean_       = param.mean_;
  sigma_      = param.sigma_;
  stat_mean_  = param.stat_mean_;
  stat_sigma_ = param.stat_sigma_;
  return *this;
}

/* resize the set of parameter */
void ModelParameters<Clust::Gaussian_sjk_>::resize(Range const& range)
{
  for (int k = mean_.begin(); k< mean_.end(); ++k)
  {
    mean_[k].resize(range) = 0.;
    sigma_[k].resize(range) = 1.;
    stat_mean_[k].resize(range);
    stat_sigma_[k].resize(range);
  }
}

/* update statistics of the parameters. */
void ModelParameters<Clust::Gaussian_sjk_>::updateStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    stat_mean_[k].update(mean_[k]);
    stat_sigma_[k].update(sigma_[k]);
  }
}
/* Set the computed statistics */
void ModelParameters<Clust::Gaussian_sjk_>::setStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    mean_[k]  = stat_mean_[k].mean();
    sigma_[k] = stat_sigma_[k].mean();
    stat_mean_[k].release();
    stat_sigma_[k].release();
  }
}
/* Release the computed statistics */
void ModelParameters<Clust::Gaussian_sjk_>::releaseStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    stat_mean_[k].release();
    stat_sigma_[k].release();
  }
}

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gaussian_sj_>::ModelParameters(int nbCluster)
                 : mean_(nbCluster), sigma_(nbCluster)
                 , stat_mean_(nbCluster), stat_sigma_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gaussian_sj_>::ModelParameters( ModelParameters const& param)
               : mean_(param.mean_), sigma_(param.sigma_)
               , stat_mean_(param.stat_mean_), stat_sigma_(param.stat_sigma_)
{}
/* destructor */
ModelParameters<Clust::Gaussian_sj_>::~ModelParameters() {}
/* copy operator.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gaussian_sj_>& ModelParameters<Clust::Gaussian_sj_>::operator=( ModelParameters const& param)
{
  mean_       = param.mean_;
  sigma_      = param.sigma_;
  stat_mean_  = param.stat_mean_;
  stat_sigma_ = param.stat_sigma_;
  return *this;
}
/* resize the set of parameter */
void ModelParameters<Clust::Gaussian_sj_>::resize(Range const& range)
{
  for (int k = mean_.begin(); k< mean_.end(); ++k)
  {
    mean_[k].resize(range) = 0;
    stat_mean_[k].resize(range);
  }
  sigma_.resize(range) = 1.;
  stat_sigma_.resize(range);
}
/* update statistics of the parameters. */
void ModelParameters<Clust::Gaussian_sj_>::updateStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    stat_mean_[k].update(mean_[k]);
  }
  stat_sigma_.update(sigma_);
}
/* Set the computed statistics */
void ModelParameters<Clust::Gaussian_sj_>::setStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    mean_[k]  = stat_mean_[k].mean();
    stat_mean_[k].release();
  }
  sigma_ = stat_sigma_.mean();
  stat_sigma_.release();
}
/* Release the computed statistics */
void ModelParameters<Clust::Gaussian_sj_>::releaseStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    stat_mean_[k].release();
  }
  stat_sigma_.release();
}
/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gaussian_sk_>::ModelParameters(int nbCluster)
               : mean_(nbCluster), sigma_(nbCluster)
               , stat_mean_(nbCluster), stat_sigma_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gaussian_sk_>::ModelParameters( ModelParameters const& param)
               : mean_(param.mean_), sigma_(param.sigma_)
               , stat_mean_(param.stat_mean_), stat_sigma_(param.stat_sigma_)
{}
/* destructor */
ModelParameters<Clust::Gaussian_sk_>::~ModelParameters() {}
/* copy operator.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gaussian_sk_>& ModelParameters<Clust::Gaussian_sk_>::operator=( ModelParameters const& param)
{
  mean_       = param.mean_;
  sigma_      = param.sigma_;
  stat_mean_  = param.stat_mean_;
  stat_sigma_ = param.stat_sigma_;
  return *this;
}
/* resize the set of parameter */
void ModelParameters<Clust::Gaussian_sk_>::resize(Range const& range)
{
  for (int k = mean_.begin(); k< mean_.end(); ++k)
  {
    mean_[k].resize(range) =0;
    sigma_[k] = 1.;
    stat_mean_[k].resize(range);
    stat_sigma_[k].release();
  }
}
/* update statistics of the parameters. */
void ModelParameters<Clust::Gaussian_sk_>::updateStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    stat_mean_[k].update(mean_[k]);
    stat_sigma_[k].update(sigma_[k]);
  }
}
/* Set the computed statistics */
void ModelParameters<Clust::Gaussian_sk_>::setStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    mean_[k]  = stat_mean_[k].mean();
    sigma_[k] = stat_sigma_[k].mean();
    stat_mean_[k].release();
    stat_sigma_[k].release();
  }
}
/* Release the computed statistics */
void ModelParameters<Clust::Gaussian_sk_>::releaseStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    stat_mean_[k].release();
    stat_sigma_[k].release();
  }
}
/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gaussian_s_>::ModelParameters(int nbCluster)
               : mean_(nbCluster), sigma_(1.)
               , stat_mean_(nbCluster), stat_sigma_()
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gaussian_s_>::ModelParameters( ModelParameters const& param)
               : mean_(param.mean_), sigma_(param.sigma_)
               , stat_mean_(param.stat_mean_), stat_sigma_(param.stat_sigma_)
{}
/* destructor */
ModelParameters<Clust::Gaussian_s_>::~ModelParameters() {}
/* copy operator.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gaussian_s_>& ModelParameters<Clust::Gaussian_s_>::operator=( ModelParameters const& param)
{
  mean_  = param.mean_;
  sigma_ = param.sigma_;
  stat_mean_  = param.stat_mean_;
  stat_sigma_ = param.stat_sigma_;
  return *this;
}
/* resize the set of parameter */
void ModelParameters<Clust::Gaussian_s_>::resize(Range const& range)
{
  for (int k = mean_.begin(); k< mean_.end(); ++k)
  {
    mean_[k].resize(range) = 0;
    stat_mean_[k].resize(range);
  }
  sigma_ = 1.;
  stat_sigma_.release();
}
/* update statistics of the parameters. */
void ModelParameters<Clust::Gaussian_s_>::updateStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    stat_mean_[k].update(mean_[k]);
  }
  stat_sigma_.update(sigma_);
}
/* Set the computed statistics */
void ModelParameters<Clust::Gaussian_s_>::setStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    mean_[k]  = stat_mean_[k].mean();
    stat_mean_[k].release();
  }
  sigma_ = stat_sigma_.mean();
  stat_sigma_.release();
}
/* Release the computed statistics */
void ModelParameters<Clust::Gaussian_s_>::releaseStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  { stat_mean_[k].release();}
  stat_sigma_.release();
}

} // namespace STK

