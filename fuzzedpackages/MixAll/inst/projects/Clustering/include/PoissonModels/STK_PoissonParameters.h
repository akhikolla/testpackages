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
 *  @brief In this file we define the Parameters classes for Poisson
 *  mixture models
 **/

#ifndef STK_POISSONPARAMETERS_H
#define STK_POISSONPARAMETERS_H

#include "../STK_Clust_Util.h"

#include <Arrays/include/STK_Array1D.h>
#include <Arrays/include/STK_CArrayPoint.h>
#include <Arrays/include/STK_CArrayVector.h>

#include <STatistiK/include/STK_Stat_Online.h>
#include <STatistiK/include/STK_Stat_Functors.h>

namespace STK
{
/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Poisson_ljlk model.
 */
template<>
struct ModelParameters<Clust::Poisson_ljlk_>
{
    /** intensity of the variables by class */
    CVectorX lambdak_;
    /** intensity of the variables by variables */
    CPointX lambdaj_;
    /** Array of the lambdak_ statistics */
    Array1D< Stat::Online<Real, Real> > stat_lambdak_;
    /** Array of the lambdaj_ statistics */
    Stat::Online<CVectorX, Real>  stat_lambdaj_;

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

    /** @return the intensity of the kth cluster and jth variable */
    inline Real lambda(int k, int j) const { return lambdak_[k] * lambdaj_[j];}

    /** resize the set of parameter */
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array @c params is of size @c (nbCluster, nbVariable).
     **/
    template<class Array>
    void setParameters( ExprBase<Array> const& params)
    {
      lambdak_ = Stat::meanByRow(params.asDerived());
      lambdaj_ = Stat::meanByCol(params.asDerived());
      Real cte = std::sqrt((params/(lambdak_ * lambdaj_)).mean());
      lambdak_ *= cte;
      lambdaj_ *= cte;
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Poisson_ljlk model.
 */
template<>
struct ModelParameters<Clust::Poisson_ljk_>
{
    /** intensity of the variables */
    Array1D<CPointX> lambda_;
    /** Array of the lambdak_ statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_lambda_;

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

    /** @return the intensity of the kth cluster and jth variable */
    inline Real const& lambda(int k, int j) const { return lambda_[k][j];}

    /** resize the set of parameter */
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array @c params is of size @c (nbCluster, nbVariable).
     **/
    template<class Array>
    void setParameters( ExprBase<Array> const& params)
    {
      for(int k=lambda_.begin(); k<lambda_.end(); ++k)
      { lambda_[k] = params.row(k);}
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Poisson_lk model.
 */
template<>
struct ModelParameters<Clust::Poisson_lk_>
{
    /** intensity of the variables */
    Array1D<Real> lambda_;
    /** Array of the lambdak_ statistics */
    Array1D< Stat::Online<Real, Real> > stat_lambda_;

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

    /** @return the intensity of the kth cluster and jth variable */
    inline Real const& lambda(int k, int j) const { return lambda_[k];}

    /** resize the set of parameter */
    void resize(Range const& range);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array @c params is of size @c (nbCluster, nbVariable).
     **/
    template<class Array>
    void setParameters( ExprBase<Array> const& params)
    {
      for(int k=lambda_.begin(); k<lambda_.end(); ++k)
      { lambda_[k] = params.row(k).mean();}
    }
};


} // namespace STK

#endif /* STK_POISSONPARAMETERS_H */
