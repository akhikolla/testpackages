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

/** @file  STK_GammaParameters.h
 *  @brief In this file we define the Parameters classes for gamma
 *  mixture models
 **/

#ifndef STK_GAMMAPARAMETERS_H
#define STK_GAMMAPARAMETERS_H

#include "../STK_Clust_Util.h"

#include <Arrays/include/STK_Array1D.h>
#include <Arrays/include/STK_CArrayPoint.h>

#include <STatistiK/include/STK_Stat_Online.h>

namespace STK
{
/** base class of the Gamma models Parameter Handler */
struct ParametersGammaBase
{
  /** mean */
  Array1D<CPointX> mean_;
  /** log-means */
  Array1D<CPointX> meanLog_;
  /** variance_ */
  Array1D<CPointX> variance_;

  /** default constructor */
  ParametersGammaBase( int nbCluster);
  /** copy constructor */
  ParametersGammaBase( ParametersGammaBase const& model);
  /** destructor */
  ~ParametersGammaBase();
  /** copy operator */
  ParametersGammaBase& operator=( ParametersGammaBase const& other);
  /** @param range the range of the variables */
  void resize(Range const& range);
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Gamma_a_bjk_ mixture model.
 */
template<>
struct ModelParameters<Clust::Gamma_a_bjk_>: ParametersGammaBase
{
    /** shape of the variables */
    Real shape_;
    /** scales of the variables */
    Array1D<CPointX>  scale_;
    /** shape statistics */
    Stat::Online<Real, Real> stat_shape_;
    /** Array of the scale statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_scale_;

   /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param);
    /** destructor */
    ~ModelParameters();

    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_;}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_[k][j];}

    /** resize and initialize the set of parameter.
     *  @param range range of the variables in the data set
     **/
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters( ExprBase<Array> const& params)
    {
      shape_ = 0.;
      for(int k=scale_.begin(), kp= params.beginRows(); k<scale_.end(); ++k, kp+=2)
      {
        shape_ += params.row(kp).mean();
        scale_[k] = params.row(kp+1);
      }
      shape_ /= scale_.size();
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Gamma_a_bk_ mixture model.
 */
template<>
struct ModelParameters<Clust::Gamma_a_bk_>: public ParametersGammaBase
{
    /** shape of the variables */
    Real shape_;
    /** scales of the variables */
    Array1D<Real>  scale_;
    /** shape statistics */
    Stat::Online<Real, Real> stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_scale_;

    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param);
    /** destructor */
    ~ModelParameters();

    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_;}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_[k];}

    /** resize and initialize the set of parameter.
     *  @param range range of the variables in the data set
     **/
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters(ExprBase<Array> const& params)
    {
      shape_ = 0.;
      for(int k=scale_.begin(), kp= params.beginRows(); k<scale_.end(); ++k, kp+=2)
      {
        shape_ += params.row(kp).mean();
        scale_[k] = params.row(kp+1).mean();
      }
      shape_ /= scale_.size();
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Gamma_aj_bjk_ mixture model.
 */
template<>
struct ModelParameters<Clust::Gamma_aj_bjk_>: public ParametersGammaBase
{
    /** shapes of the variables */
    CPointX shape_;
    /** scales of the variables */
    Array1D<CPointX>  scale_;
    /** shape statistics */
    Stat::Online<CPointX, Real> stat_shape_;
    /** Array of the scale statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_scale_;

    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param);
    /** destructor */
    ~ModelParameters();

    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_[k];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_[k][j];}

    /** resize and initialize the set of parameter.
     *  @param range range of the variables in the data set
     **/
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters(ExprBase<Array> const& params)
    {
      shape_ = 0.;
      for(int k=scale_.begin(), kp= params.beginRows(); k<scale_.end(); ++k, kp+=2)
      {
        shape_ += params.row(kp);
        scale_[k] = params.row(kp+1);
      }
      shape_ /= scale_.size();
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Gamma_aj_bk_ mixture model.
 */
template<>
struct ModelParameters<Clust::Gamma_aj_bk_>: public ParametersGammaBase
{
    /** shapes of the variables */
    CPointX shape_;
    /** scales of the variables */
    Array1D<Real>  scale_;
    /** shape statistics */
    Stat::Online<CPointX, Real> stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_scale_;

    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param);
    /** destructor */
    ~ModelParameters();

    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_[j];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_[k];}

    /** resize and initialize the set of parameter.
     *  @param range range of the variables in the data set
     **/
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters(ExprBase<Array> const& params)
    {
      shape_ = 0.;
      for(int k=scale_.begin(), kp= params.beginRows(); k<scale_.end(); ++k, kp+=2)
      {
        shape_ += params.row(kp);
        scale_[k] = params.row(kp+1).mean();
      }
      shape_ /= scale_.size();
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Gamma_ajk_b_ mixture model.
 */
template<>
struct ModelParameters<Clust::Gamma_ajk_b_>: public ParametersGammaBase
{
    /** shapes of the variables */
    Array1D<CPointX> shape_;
    /** scales of the variables */
    Real  scale_;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Stat::Online<Real, Real> stat_scale_;

    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param);
    /** destructor */
    ~ModelParameters();

    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_[k][j];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_;}

    /** resize and initialize the set of parameter.
     *  @param range range of the variables in the data set
     **/
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters(ExprBase<Array> const& params)
    {
      scale_ = 0.;
      for(int k=shape_.begin(), kp= params.beginRows(); k<shape_.end(); ++k, kp+=2)
      {
        shape_[k] = params.row(kp);
        scale_ = params.row(kp+1).mean();
      }
      scale_ /= shape_.size();
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Gamma_ajk_bj_ mixture model.
 */
template<>
struct ModelParameters<Clust::Gamma_ajk_bj_>: public ParametersGammaBase
{
    /** shapes of the variables */
    Array1D<CPointX> shape_;
    /** scales of the variables */
    CPointX  scale_;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Stat::Online<CPointX, Real> stat_scale_;

    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param);
    /** destructor */
    ~ModelParameters();

    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_[k][j];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_[j];}

    /** resize and initialize the set of parameter.
     *  @param range range of the variables in the data set
     **/
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters(ExprBase<Array> const& params)
    {
      scale_ = 0.;
      for(int k=shape_.begin(), kp= params.beginRows(); k<shape_.end(); ++k, kp+=2)
      {
        shape_[k] = params.row(kp);
        scale_ += params.row(kp+1);
      }
      scale_ /= shape_.size();
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Gamma_ajk_bjk_ mixture model.
 */
template<>
struct ModelParameters<Clust::Gamma_ajk_bjk_>: public ParametersGammaBase
{
    /** shapes of the variables */
    Array1D<CPointX> shape_;
    /** scales of the variables */
    Array1D<CPointX> scale_;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_scale_;

    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param);
    /** destructor */
    ~ModelParameters();

    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_[k][j];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_[k][j];}

    /** resize and initialize the set of parameter.
     *  @param range range of the variables in the data set
     **/
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters(ExprBase<Array> const& params)
    {
      for(int k=shape_.begin(), kp= params.beginRows(); k<shape_.end(); ++k, kp+=2)
      {
        shape_[k] = params.row(kp);
        scale_[k] = params.row(kp+1);
      }
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Gamma_ajk_bk_ mixture model.
 */
template<>
struct ModelParameters<Clust::Gamma_ajk_bk_>: public ParametersGammaBase
{
    /** shapes of the variables */
    Array1D<CPointX> shape_;
    /** scales of the variables */
    Array1D<Real> scale_;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_scale_;

    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param);
    /** destructor */
    ~ModelParameters();

    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_[k][j];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_[k];}

    /** resize and initialize the set of parameter.
     *  @param range range of the variables in the data set
     **/
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters(ExprBase<Array> const& params)
    {
      for(int k=shape_.begin(), kp= params.beginRows(); k<shape_.end(); ++k, kp+=2)
      {
        shape_[k] = params.row(kp);
        scale_[k] = params.row(kp+1).mean();
      }
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Gamma_ak_b_ mixture model.
 */
template<>
struct ModelParameters<Clust::Gamma_ak_b_>: public ParametersGammaBase
{
    /** shapes of the variables */
    Array1D<Real> shape_;
    /** scales of the variables */
    Real  scale_;
    /** Array of the mean statistics */
    Array1D< Stat::Online<Real, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Stat::Online<Real, Real> stat_scale_;

    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param);
    /** destructor */
    ~ModelParameters();

    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_[k];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_;}

    /** resize and initialize the set of parameter.
     *  @param range range of the variables in the data set
     **/
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters(ExprBase<Array> const& params)
    {
      scale_ = 0.;
      for(int k=shape_.begin(), kp= params.beginRows(); k<shape_.end(); ++k, kp+=2)
      {
        shape_[k] = params.row(kp).mean();
        scale_ += params.row(kp+1).mean();
      }
      scale_ /= shape_.size();
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Gamma_ak_bj_ mixture model.
 */
template<>
struct ModelParameters<Clust::Gamma_ak_bj_>: public ParametersGammaBase
{
    /** shapes of the variables */
    Array1D<Real> shape_;
    /** scales of the variables */
    CPointX  scale_;
    /** Array of the mean statistics */
    Array1D< Stat::Online<Real, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Stat::Online<CPointX, Real> stat_scale_;

    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param);
    /** destructor */
    ~ModelParameters();

    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_[k];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_[j];}

    /** resize and initialize the set of parameter.
     *  @param range range of the variables in the data set
     **/
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters(ExprBase<Array> const& params)
    {
      scale_ = 0.;
      for(int k=shape_.begin(), kp= params.beginRows(); k<shape_.end(); ++k, kp+=2)
      {
        shape_[k] = params.row(kp).mean();
        scale_ += params.row(kp+1);
      }
      scale_ /= shape_.size();
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Gamma_ak_bjk_ mixture model.
 */
template<>
struct ModelParameters<Clust::Gamma_ak_bjk_>: public ParametersGammaBase
{
    /** shapes of the variables */
    Array1D<Real> shape_;
    /** scales of the variables */
    Array1D<CPointX> scale_;
    /** Array of the mean statistics */
    Array1D< Stat::Online<Real, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_scale_;

    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param);
    /** destructor */
    ~ModelParameters();

    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_[k];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_[k][j];}

    /** resize and initialize the set of parameter.
     *  @param range range of the variables in the data set
     **/
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters(ExprBase<Array> const& params)
    {
      for(int k=shape_.begin(), kp= params.beginRows(); k<shape_.end(); ++k, kp+=2)
      {
        shape_[k] = params.row(kp).mean();
        scale_[k] = params.row(kp+1);
      }
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Gamma_ak_bk_ mixture model.
 */
template<>
struct ModelParameters<Clust::Gamma_ak_bk_>: public ParametersGammaBase
{
    /** shapes of the variables */
    Array1D<Real> shape_;
    /** scales of the variables */
    Array1D<Real> scale_;
    /** Array of the mean statistics */
    Array1D< Stat::Online<Real, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_scale_;

    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param);
    /** destructor */
    ~ModelParameters();

    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_[k];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_[k];}

    /** resize and initialize the set of parameter.
     *  @param range range of the variables in the data set
     **/
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters(ExprBase<Array> const& params)
    {
      for(int k=shape_.begin(), kp= params.beginRows(); k<shape_.end(); ++k, kp+=2)
      {
        shape_[k] = params.row(kp).mean();
        scale_[k] = params.row(kp+1).mean();
      }
    }
};

} // namespace STK

#endif /* STK_GAMMAPARAMETERS_H */
