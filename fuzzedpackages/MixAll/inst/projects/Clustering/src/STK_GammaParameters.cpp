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

/** @file  STK_GammaParameters.cpp
 *  @brief In this file we implement the Parameters classes for gamma
 *  mixture models
 **/

#include "../include/GammaModels/STK_GammaParameters.h"

namespace STK
{
/* copy operator */
ParametersGammaBase& ParametersGammaBase::operator=( ParametersGammaBase const& other)
{
  mean_     = other.mean_;
  meanLog_  = other.meanLog_;
  variance_ = other.variance_;
  return *this;
}

/* default constructor */
ParametersGammaBase::ParametersGammaBase( int nbCluster)
                   : mean_(nbCluster), meanLog_(nbCluster)
                   , variance_(nbCluster)
{}
/* copy constructor */
ParametersGammaBase::ParametersGammaBase( ParametersGammaBase const& model)
                   : mean_(model.mean_)
                   , meanLog_(model.meanLog_)
                   , variance_(model.variance_) {}
/* destructor */
ParametersGammaBase::~ParametersGammaBase() {}
/* @param range the range of the variables */
void ParametersGammaBase::resize(Range const& range)
{
  for (int k=mean_.begin(); k < mean_.end(); ++k)
  {
    mean_[k].resize(range) = 1.;
    meanLog_[k].resize(range) = 0.;
    variance_[k].resize(range) = 1.;
  }
}

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gamma_a_bjk_>::ModelParameters(int nbCluster)
              : ParametersGammaBase(nbCluster)
              , shape_(0.), scale_(nbCluster)
              , stat_shape_(), stat_scale_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gamma_a_bjk_>::ModelParameters( ModelParameters const& param)
               : ParametersGammaBase(param)
               , shape_(param.shape_), scale_(param.scale_)
               , stat_shape_(param.stat_shape_), stat_scale_(param.stat_scale_)
{}
/* destructor */
ModelParameters<Clust::Gamma_a_bjk_>::~ModelParameters() {}

/* resize the set of parameter */
void ModelParameters<Clust::Gamma_a_bjk_>::resize(Range const& range)
{
  ParametersGammaBase::resize(range);
  shape_ = 1.;
  stat_shape_.release();
  for (int k = scale_.begin(); k< scale_.end(); ++k)
  {
    scale_[k].resize(range) = 1.;
    stat_scale_[k].resize(range);
  }
}

/* update statistics of the parameters. */
void ModelParameters<Clust::Gamma_a_bjk_>::updateStatistics()
{
  for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
  { stat_scale_[k].update(scale_[k]);}
  stat_shape_.update(shape_);
}
/* Set the computed statistics */
void ModelParameters<Clust::Gamma_a_bjk_>::setStatistics()
{
  for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
  {
    scale_[k] = stat_scale_[k].mean();
    stat_scale_[k].release();
  }
  shape_ = stat_shape_.mean();
  stat_shape_.release();
}
/* Release the computed statistics */
void ModelParameters<Clust::Gamma_a_bjk_>::releaseStatistics()
{
  for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
  { stat_scale_[k].release();}
  stat_shape_.release();
}

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gamma_a_bk_>::ModelParameters(int nbCluster)
              : ParametersGammaBase(nbCluster)
              , shape_(0.), scale_(nbCluster)
              , stat_shape_(), stat_scale_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gamma_a_bk_>::ModelParameters( ModelParameters const& param)
               : ParametersGammaBase(param)
               , shape_(param.shape_), scale_(param.scale_)
               , stat_shape_(param.stat_shape_), stat_scale_(param.stat_scale_)
{}
/* destructor */
ModelParameters<Clust::Gamma_a_bk_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Gamma_a_bk_>::resize(Range const& range)
{
  ParametersGammaBase::resize(range);
  shape_ = 1.;
  stat_shape_.release();
  for (int k = scale_.begin(); k< scale_.end(); ++k)
  { scale_[k] = 1.;}
}

/* update statistics of the parameters. */
void ModelParameters<Clust::Gamma_a_bk_>::updateStatistics()
{
  for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
  { stat_scale_[k].update(scale_[k]);}
  stat_shape_.update(shape_);
}
/* Set the computed statistics */
void ModelParameters<Clust::Gamma_a_bk_>::setStatistics()
{
  for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
  {
    scale_[k] = stat_scale_[k].mean();
    stat_scale_[k].release();
  }
  shape_ = stat_shape_.mean();
  stat_shape_.release();
}
/* Release the computed statistics */
void ModelParameters<Clust::Gamma_a_bk_>::releaseStatistics()
{
  for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
  { stat_scale_[k].release();}
  stat_shape_.release();
}

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gamma_aj_bjk_>::ModelParameters(int nbCluster)
               : ParametersGammaBase(nbCluster)
               , shape_(), scale_(nbCluster)
               , stat_shape_(), stat_scale_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gamma_aj_bjk_>::ModelParameters( ModelParameters const& param)
               : ParametersGammaBase(param)
               , shape_(param.shape_), scale_(param.scale_)
               , stat_shape_(param.stat_shape_), stat_scale_(param.stat_scale_)
{}
/* destructor */
ModelParameters<Clust::Gamma_aj_bjk_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Gamma_aj_bjk_>::resize(Range const& range)
{
  ParametersGammaBase::resize(range);
  shape_.resize(range) = 1.;
  stat_shape_.resize(range);
  for (int k = scale_.begin(); k< scale_.end(); ++k)
  {
    scale_[k].resize(range) = 1.;
    stat_scale_[k].resize(range);
  }
}

/* update statistics of the parameters. */
void ModelParameters<Clust::Gamma_aj_bjk_>::updateStatistics()
{
  stat_shape_.release();
  for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
  { stat_scale_[k].update(scale_[k]);}
}
/* Set the computed statistics */
void ModelParameters<Clust::Gamma_aj_bjk_>::setStatistics()
{
  shape_ = stat_shape_.mean();
  stat_shape_.release();
  for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
  {
    scale_[k] = stat_scale_[k].mean();
    stat_scale_[k].release();
  }
}
/* Release the computed statistics */
void ModelParameters<Clust::Gamma_aj_bjk_>::releaseStatistics()
{
  stat_shape_.release();
  for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
  { stat_scale_[k].release();}
}

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gamma_aj_bk_>::ModelParameters(int nbCluster)
               : ParametersGammaBase(nbCluster)
               , shape_(), scale_(nbCluster)
               , stat_shape_(nbCluster), stat_scale_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gamma_aj_bk_>::ModelParameters( ModelParameters const& param)
               : ParametersGammaBase(param)
               , shape_(param.shape_), scale_(param.scale_)
               , stat_shape_(param.stat_shape_), stat_scale_(param.stat_scale_)
{}
/* destructor */
ModelParameters<Clust::Gamma_aj_bk_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Gamma_aj_bk_>::resize(Range const& range)
{
  ParametersGammaBase::resize(range);
  shape_.resize(range) = 1.;
  stat_shape_.resize(range);
  for (int k = scale_.begin(); k< scale_.end(); ++k)
  {
    scale_[k] = 1.;
    stat_scale_[k].release();
  }
}

/* update statistics of the parameters. */
void ModelParameters<Clust::Gamma_aj_bk_>::updateStatistics()
{
  for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
  { stat_scale_[k].update(scale_[k]);}
  stat_shape_.update(shape_);
}
/* Set the computed statistics */
void ModelParameters<Clust::Gamma_aj_bk_>::setStatistics()
{
  for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
  {
    scale_[k] = stat_scale_[k].mean();
    stat_scale_[k].release();
  }
  shape_ = stat_shape_.mean();
  stat_shape_.release();
}
/* Release the computed statistics */
void ModelParameters<Clust::Gamma_aj_bk_>::releaseStatistics()
{
  for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
  { stat_scale_[k].release();}
  stat_shape_.release();
}

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gamma_ajk_b_>::ModelParameters(int nbCluster)
               : ParametersGammaBase(nbCluster)
               , shape_(nbCluster), scale_()
               , stat_shape_(nbCluster), stat_scale_()
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gamma_ajk_b_>::ModelParameters( ModelParameters const& param)
               : ParametersGammaBase(param)
               , shape_(param.shape_), scale_(param.scale_)
               , stat_shape_(param.stat_shape_), stat_scale_(param.stat_scale_)
{}
/* destructor */
ModelParameters<Clust::Gamma_ajk_b_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Gamma_ajk_b_>::resize(Range const& range)
{
  ParametersGammaBase::resize(range);
  for (int k = shape_.begin(); k< shape_.end(); ++k)
  {
    shape_[k].resize(range) = 1.;
    stat_shape_[k].resize(range);
  }
  scale_ = 1.;
  stat_scale_.release();
}
/* update statistics of the parameters. */
void ModelParameters<Clust::Gamma_ajk_b_>::updateStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  { stat_shape_[k].update(shape_[k]);}
  stat_scale_.update(scale_);
}
/* Set the computed statistics */
void ModelParameters<Clust::Gamma_ajk_b_>::setStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  {
    shape_[k] = stat_shape_[k].mean();
    stat_shape_[k].release();
  }
  scale_ = stat_scale_.mean();
  stat_scale_.release();
}
/* Release the computed statistics */
void ModelParameters<Clust::Gamma_ajk_b_>::releaseStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  { stat_shape_[k].release();}
  stat_scale_.release();
}

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gamma_ajk_bj_>::ModelParameters(int nbCluster)
               : ParametersGammaBase(nbCluster)
               , shape_(nbCluster), scale_()
               , stat_shape_(nbCluster), stat_scale_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gamma_ajk_bj_>::ModelParameters( ModelParameters const& param)
               : ParametersGammaBase(param)
               , shape_(param.shape_), scale_(param.scale_)
               , stat_shape_(param.stat_shape_), stat_scale_(param.stat_scale_)
{}
/* destructor */
ModelParameters<Clust::Gamma_ajk_bj_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Gamma_ajk_bj_>::resize(Range const& range)
{
  ParametersGammaBase::resize(range);
  scale_.resize(range) = 1.;
  stat_scale_.resize(range);
  for (int k = shape_.begin(); k< shape_.end(); ++k)
  {
    shape_[k].resize(range) = 1.;
    stat_shape_[k].resize(range);
  }
}
/* update statistics of the parameters. */
void ModelParameters<Clust::Gamma_ajk_bj_>::updateStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  { stat_shape_[k].update(shape_[k]);}
  stat_scale_.update(scale_);
}
/* Set the computed statistics */
void ModelParameters<Clust::Gamma_ajk_bj_>::setStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  {
    shape_[k] = stat_shape_[k].mean();
    stat_shape_[k].release();
  }
  scale_ = stat_scale_.mean();
  stat_scale_.release();
}
/* Release the computed statistics */
void ModelParameters<Clust::Gamma_ajk_bj_>::releaseStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  { stat_shape_[k].release();}
  stat_scale_.release();
}


/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gamma_ajk_bjk_>::ModelParameters(int nbCluster)
               : ParametersGammaBase(nbCluster)
               , shape_(nbCluster), scale_(nbCluster)
               , stat_shape_(nbCluster), stat_scale_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gamma_ajk_bjk_>::ModelParameters( ModelParameters const& param)
               : ParametersGammaBase(param)
               , shape_(param.shape_), scale_(param.scale_)
               , stat_shape_(param.stat_shape_), stat_scale_(param.stat_scale_)
{}
/* destructor */
ModelParameters<Clust::Gamma_ajk_bjk_>::~ModelParameters() {}

/* resize the set of parameter */
void ModelParameters<Clust::Gamma_ajk_bjk_>::resize(Range const& range)
{
  ParametersGammaBase::resize(range);
  for (int k = shape_.begin(); k< shape_.end(); ++k)
  {
    shape_[k].resize(range) = 1.;
    stat_shape_[k].resize(range);
    scale_[k].resize(range) = 1.;
    stat_scale_[k].resize(range);
  }
}

/* update statistics of the parameters. */
void ModelParameters<Clust::Gamma_ajk_bjk_>::updateStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  { stat_shape_[k].update(shape_[k]);
    stat_scale_[k].update(scale_[k]);
  }
}
/* Set the computed statistics */
void ModelParameters<Clust::Gamma_ajk_bjk_>::setStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  {
    shape_[k] = stat_shape_[k].mean();
    stat_shape_[k].release();
    scale_[k] = stat_scale_[k].mean();
    stat_scale_[k].release();
  }
}
/* Release the computed statistics */
void ModelParameters<Clust::Gamma_ajk_bjk_>::releaseStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  {
    stat_shape_[k].release();
    stat_scale_[k].release();
  }
}

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gamma_ajk_bk_>::ModelParameters(int nbCluster)
               : ParametersGammaBase(nbCluster)
               , shape_(nbCluster), scale_(nbCluster)
               , stat_shape_(nbCluster), stat_scale_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gamma_ajk_bk_>::ModelParameters( ModelParameters const& param)
               : ParametersGammaBase(param)
               , shape_(param.shape_), scale_(param.scale_)
               , stat_shape_(param.stat_shape_), stat_scale_(param.stat_scale_)
{}
/* destructor */
ModelParameters<Clust::Gamma_ajk_bk_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Gamma_ajk_bk_>::resize(Range const& range)
{
  ParametersGammaBase::resize(range);
  for (int k = shape_.begin(); k< shape_.end(); ++k)
  {
    shape_[k].resize(range) = 1.;
    stat_shape_[k].resize(range);
    scale_[k] = 1.;
    stat_scale_[k].release();
  }
}

/* update statistics of the parameters. */
void ModelParameters<Clust::Gamma_ajk_bk_>::updateStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  {
    stat_shape_[k].update(shape_[k]);
    stat_scale_[k].update(scale_[k]);
  }
}
/* Set the computed statistics */
void ModelParameters<Clust::Gamma_ajk_bk_>::setStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  {
    shape_[k] = stat_shape_[k].mean();
    stat_shape_[k].release();
    scale_[k] = stat_scale_[k].mean();
    stat_scale_[k].release();
  }
}
/* Release the computed statistics */
void ModelParameters<Clust::Gamma_ajk_bk_>::releaseStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  {
    stat_shape_[k].release();
    stat_scale_[k].release();
  }
}

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gamma_ak_b_>::ModelParameters(int nbCluster)
               : ParametersGammaBase(nbCluster)
               , shape_(nbCluster), scale_()
               , stat_shape_(nbCluster), stat_scale_()
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gamma_ak_b_>::ModelParameters( ModelParameters const& param)
               : ParametersGammaBase(param)
               , shape_(param.shape_), scale_(param.scale_)
               , stat_shape_(param.stat_shape_), stat_scale_(param.stat_scale_)
{}
/* destructor */
ModelParameters<Clust::Gamma_ak_b_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Gamma_ak_b_>::resize(Range const& range)
{
  ParametersGammaBase::resize(range);
  scale_ = 1.;
  stat_scale_.release();
  for (int k = shape_.begin(); k< shape_.end(); ++k)
  {
    shape_[k] = 1.;
    stat_shape_[k].release();
  }
}

/* update statistics of the parameters. */
void ModelParameters<Clust::Gamma_ak_b_>::updateStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  { stat_shape_[k].update(shape_[k]);}
  stat_scale_.update(scale_);
}
/* Set the computed statistics */
void ModelParameters<Clust::Gamma_ak_b_>::setStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  {
    shape_[k] = stat_shape_[k].mean();
    stat_shape_[k].release();
  }
  scale_ = stat_scale_.mean();
  stat_scale_.release();
}
/* Release the computed statistics */
void ModelParameters<Clust::Gamma_ak_b_>::releaseStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  { stat_shape_[k].release();}
  stat_scale_.release();
}

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gamma_ak_bj_>::ModelParameters(int nbCluster)
               : ParametersGammaBase(nbCluster)
               , shape_(nbCluster), scale_()
               , stat_shape_(nbCluster), stat_scale_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gamma_ak_bj_>::ModelParameters( ModelParameters const& param)
               : ParametersGammaBase(param)
               , shape_(param.shape_), scale_(param.scale_)
               , stat_shape_(param.stat_shape_), stat_scale_(param.stat_scale_)
{}
/* destructor */
ModelParameters<Clust::Gamma_ak_bj_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Gamma_ak_bj_>::resize(Range const& range)
{
  ParametersGammaBase::resize(range);
  for (int k = shape_.begin(); k< shape_.end(); ++k)
  {
    shape_[k] = 1.;
    stat_shape_[k].release();
  }
  scale_.resize(range) = 1.;
  stat_scale_.resize(range);
}

/* update statistics of the parameters. */
void ModelParameters<Clust::Gamma_ak_bj_>::updateStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  { stat_shape_[k].update(shape_[k]);}
  stat_scale_.update(scale_);
}
/* Set the computed statistics */
void ModelParameters<Clust::Gamma_ak_bj_>::setStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  {
    shape_[k] = stat_shape_[k].mean();
    stat_shape_[k].release();
  }
  scale_ = stat_scale_.mean();
  stat_scale_.release();
}
/* Release the computed statistics */
void ModelParameters<Clust::Gamma_ak_bj_>::releaseStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  { stat_shape_[k].release();}
  stat_scale_.release();
}

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gamma_ak_bjk_>::ModelParameters(int nbCluster)
               : ParametersGammaBase(nbCluster)
               , shape_(nbCluster), scale_(nbCluster)
               , stat_shape_(nbCluster), stat_scale_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gamma_ak_bjk_>::ModelParameters( ModelParameters const& param)
               : ParametersGammaBase(param)
               , shape_(param.shape_), scale_(param.scale_)
               , stat_shape_(param.stat_shape_), stat_scale_(param.stat_scale_)
{}
/* destructor */
ModelParameters<Clust::Gamma_ak_bjk_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Gamma_ak_bjk_>::resize(Range const& range)
{
  ParametersGammaBase::resize(range);
  for (int k = shape_.begin(); k< shape_.end(); ++k)
  {
    shape_[k] = 1.;
    stat_shape_[k].release();
    scale_[k].resize(range) = 1.;
    stat_scale_[k].resize(range);
  }
}
/* update statistics of the parameters. */
void ModelParameters<Clust::Gamma_ak_bjk_>::updateStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  { stat_shape_[k].update(shape_[k]);
    stat_scale_[k].update(scale_[k]);
  }
}
/* Set the computed statistics */
void ModelParameters<Clust::Gamma_ak_bjk_>::setStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  {
    shape_[k] = stat_shape_[k].mean();
    stat_shape_[k].release();
    scale_[k] = stat_scale_[k].mean();
    stat_scale_[k].release();
  }
}
/* Release the computed statistics */
void ModelParameters<Clust::Gamma_ak_bjk_>::releaseStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  {
    stat_shape_[k].release();
    stat_scale_[k].release();
  }
}


/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Gamma_ak_bk_>::ModelParameters(int nbCluster)
               : ParametersGammaBase(nbCluster)
               , shape_(nbCluster), scale_(nbCluster)
               , stat_shape_(nbCluster), stat_scale_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Gamma_ak_bk_>::ModelParameters( ModelParameters const& param)
               : ParametersGammaBase(param)
               , shape_(param.shape_), scale_(param.scale_)
               , stat_shape_(param.stat_shape_), stat_scale_(param.stat_scale_)
{}
/* destructor */
ModelParameters<Clust::Gamma_ak_bk_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Gamma_ak_bk_>::resize(Range const& range)
{
  ParametersGammaBase::resize(range);
  for (int k = shape_.begin(); k< shape_.end(); ++k)
  {
    shape_[k] = 1.;
    stat_shape_[k].release();
    scale_[k] = 1.;
    stat_scale_[k].release();
  }
}
/* update statistics of the parameters. */
void ModelParameters<Clust::Gamma_ak_bk_>::updateStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  { stat_shape_[k].update(shape_[k]);
    stat_scale_[k].update(scale_[k]);
  }
}
/* Set the computed statistics */
void ModelParameters<Clust::Gamma_ak_bk_>::setStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  {
    shape_[k] = stat_shape_[k].mean();
    stat_shape_[k].release();
    scale_[k] = stat_scale_[k].mean();
    stat_scale_[k].release();
  }
}
/* Release the computed statistics */
void ModelParameters<Clust::Gamma_ak_bk_>::releaseStatistics()
{
  for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
  {
    stat_shape_[k].release();
    stat_scale_[k].release();
  }
}

} // namespace STK

