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

/** @file  STK_DiagGaussianParameters.h
 *  @brief In this file we define the Parameters classes for Diagonal Gaussian
 *  mixture models
 **/

#ifndef STK_DIAGGAUSSIANPARAMETERS_H
#define STK_DIAGGAUSSIANPARAMETERS_H

#include "../STK_Clust_Util.h"

#include <Arrays/include/STK_Array1D.h>
#include <Arrays/include/STK_CArray.h> // for Gaussian_sksj model
#include <Arrays/include/STK_CArrayPoint.h>
#include <Arrays/include/STK_CArrayVector.h>

#include <STatistiK/include/STK_Stat_Online.h>
#include <STatistiK/include/STK_Stat_Functors.h>

namespace STK
{

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a STK::DiagGaussian_sjk model.
 */
template<>
struct ModelParameters<Clust::Gaussian_sjk_>
{
  /** array of size nbCluster with the parameters mean of the variables */
  Array1D<CPointX> mean_;
  /** standard deviation of the variables */
  Array1D<CPointX> sigma_;
  /** Array of the mean statistics */
  Array1D< Stat::Online<CPointX, Real> > stat_mean_;
  /** Array of the standard deviation statistics */
  Array1D< Stat::Online<CPointX, Real> > stat_sigma_;

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
  /** copy operator.
   *  @param param the parameters to copy.
   **/
  ModelParameters& operator=( ModelParameters const& param);

  /** @return the mean of the kth cluster and jth variable */
  inline Real const& mean(int k, int j) const { return mean_[k][j];}
  /** @return the standard deviation of the kth cluster and jth variable */
  inline Real const& sigma(int k, int j) const { return sigma_[k][j];}

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
   *  It is assumed that the array params store for each class the mean and
   *  sd parameters on two consecutive rows.
   *  The number of column of params is the number of variables.
   **/
  template<class Array>
  void setParameters( ExprBase<Array> const& params)
  {
    for(int k=mean_.begin(), kp= params.beginRows(); k<mean_.end(); ++k, kp+=2)
    {
      mean_[k] = params.row(kp);
      sigma_[k] = params.row(kp+1);
    }
  }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a STK::DiagGaussian_sk model.
 */
template<>
struct ModelParameters<Clust::Gaussian_sk_>
{
  /** array of size nbCluster with the parameters mean of the variables */
  Array1D<CPointX> mean_;
  /** standard deviation of the variables */
  Array1D<Real> sigma_;
  /** Array of the mean statistics */
  Array1D< Stat::Online<CPointX, Real> > stat_mean_;
  /** Array of the standard deviation statistics */
  Array1D< Stat::Online<Real, Real> > stat_sigma_;

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
  /** copy operator.
   *  @param param the parameters to copy.
   **/
  ModelParameters& operator=( ModelParameters const& param);

  /** @return the mean of the kth cluster and jth variable */
  inline Real const& mean(int k, int j) const { return mean_[k][j];}
  /** @return the standard deviation of the kth cluster and jth variable */
  inline Real const& sigma(int k, int j) const { return sigma_[k];}

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
   *  It is assumed that the array params store for each class the means and
   *  sd parameters on two consecutive rows.
   *  The number of column of params is the number of variables.
   **/
  template<class Array>
  void setParameters( ExprBase<Array> const& params)
  {
    for(int k=mean_.begin(), kp= params.beginRows(); k<mean_.end(); ++k, kp+=2)
    {
      mean_[k] = params.row(kp);
      sigma_[k] = params.row(kp+1).mean();
    }
  }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a STK::DiagGaussian_sj model.
 */
template<>
struct ModelParameters<Clust::Gaussian_sj_>
{
  /** array of size nbCluster with the parameters mean of the variables */
  Array1D<CPointX> mean_;
  /** standard deviation of the variables */
  CPointX sigma_;
  /** Array of the mean statistics */
  Array1D< Stat::Online<CPointX, Real> > stat_mean_;
  /** Array of the standard deviation statistics */
  Stat::Online<CPointX, Real> stat_sigma_;

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
  /** copy operator.
   *  @param param the parameters to copy.
   **/
  ModelParameters& operator=( ModelParameters const& param);

  /** @return the mean of the kth cluster and jth variable */
  inline Real const& mean(int k, int j) const { return mean_[k][j];}
  /** @return the standard deviation of the kth cluster and jth variable */
  inline Real const& sigma(int k, int j) const { return sigma_[j];}

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
   *  It is assumed that the array params store for each class the means and
   *  sd parameters on two consecutive rows.
   *  The number of column of params is the number of variables.
   **/
  template<class Array>
  void setParameters( ExprBase<Array> const& params)
  {
    sigma_ =0.;
    for(int k=mean_.begin(), kp= params.beginRows(); k<mean_.end(); ++k, kp+=2)
    {
      mean_[k] = params.row(kp);
      sigma_ += params.row(kp+1);
    }
    sigma_ /= mean_.size();
  }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a STK::DiagGaussian_s model.
 */
template<>
struct ModelParameters<Clust::Gaussian_s_>
{
  /** array of size nbCluster with the parameters mean of the variables */
  Array1D<CPointX> mean_;
  /** standard deviation of the variables */
  Real sigma_;
  /** Array of the mean statistics */
  Array1D< Stat::Online<CPointX, Real> > stat_mean_;
  /** Array of the standard deviation statistics */
  Stat::Online<Real, Real> stat_sigma_;

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
  /** copy operator.
   *  @param param the parameters to copy.
   **/
  ModelParameters& operator=( ModelParameters const& param);

  /** @return the mean of the kth cluster and jth variable */
  inline Real const& mean(int k, int j) const { return mean_[k][j];}
  /** @return the standard deviation of the kth cluster and jth variable */
  inline Real const& sigma(int k, int j) const { return sigma_;}

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
   *  It is assumed that array params stores for each class the mean and
   *  sd parameters on two consecutive rows.
   *  The number of column of params is the number of variables.
   **/
  template<class Array>
  void setParameters( ExprBase<Array> const& params)
  {
    sigma_ =0.;
    for(int k=mean_.begin(), kp= params.beginRows(); k<mean_.end(); ++k, kp+=2)
    {
      mean_[k] = params.row(kp);
      sigma_ += params.row(kp+1).mean();
    }
    sigma_ /= mean_.size();
  }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Gaussian_sjsk mixture model.
 */
template<>
struct ModelParameters<Clust::Gaussian_sjsk_>
{
    /** array of size nbCluster with the parameters mean of the variables */
    Array1D<CPointX> mean_;
    /** standard deviation of the variables */
    CVectorX sigmak_;
    /** standard deviation of the variables by variables */
    CPointX sigmaj_;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_mean_;
    /** Array of the standard deviation k statistics */
    Stat::Online<Real, Real> stat_sigmak_;
    /** Array of the standard deviation j statistics */
    Stat::Online<CPointX, Real> stat_sigmaj_;

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
    inline Real const& mean(int k, int j) const { return mean_[k][j];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real sigma(int k, int j) const { return sigmak_[k] * sigmaj_[j];}

    /** resize the set of parameter */
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array @c params is of size @c (2*nbCluster, nbVariable).
     **/
    template<class Array>
    void setParameters( ExprBase<Array> const& params)
    {
      CArrayXX sigma(mean_.range(), params.cols());
      for(int k=mean_.begin(), kp= params.beginRows(); k<mean_.end(); ++k, kp+=2)
      {
          mean_[k]     = params.row(kp);
          sigma.row(k) = params.row(kp+1);
      }
      sigmak_ = Stat::meanByRow(sigma.asDerived());
      sigmaj_ = Stat::meanByCol(sigma.asDerived());
      Real cte = std::sqrt((sigma/(sigmak_ * sigmaj_)).mean());
      sigmak_ *= cte;
      sigmaj_ *= cte;
    }
};

} // namespace STK

#endif /* STK_DIAGGAUSSIANPARAMETERS_H */
