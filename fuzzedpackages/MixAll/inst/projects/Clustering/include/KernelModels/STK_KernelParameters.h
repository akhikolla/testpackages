/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2017  Serge Iovleff, Université Lille 1, Inria

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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: Sep 8, 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_KernelParameters.h
 *  @brief In this file we define the Parameters classes for kernel
 *  mixture models.
 **/


#ifndef STK_KERNELPARAMETERS_H
#define STK_KERNELPARAMETERS_H

#include "../STK_Clust_Util.h"

#include <Arrays/include/STK_Array1D.h>
#include <Arrays/include/STK_CArrayPoint.h>

#include <STatistiK/include/STK_Stat_Online.h>

namespace STK
{
/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Kmm_s model.
 */
template<>
struct ModelParameters<Clust::Kmm_s_>
{
    /** variance of the variables */
    Real sigma2_;
    /** dimension of the gaussian kernel */
    CPointX dim_;

    /** sigma2 statistics */
    Stat::Online<Real, Real> stat_sigma2_;
    /** Array of the dim statistics */
    Array1D< Stat::Online<Real, Real> > stat_dim_;

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

    /** @return the standard deviation of the kth cluster */
    inline Real const& sigma2(int k) const { return sigma2_;}
    /** @return the dimension of the kth cluster */
    inline Real const& dim(int k) const { return dim_[k];}

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the sigma2 and
     *  dim parameters on consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(ExprBase<Array> const& params)
    {
      sigma2_ = 0.;
      for(int k=dim_.begin(); k<dim_.end(); ++k)
      { sigma2_  += params(k, baseIdx)   ;
        dim_[k]  = params(k, baseIdx+1);
      }
      sigma2_ /= dim_.size();
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Kmm_sk model.
 */
template<>
struct ModelParameters<Clust::Kmm_sk_>
{
    /** variance of the variables */
    CPointX sigma2_;
    /** dimension of the gaussian kernel */
    CPointX dim_;
    /** Array of the sigma2 statistics */
    Array1D< Stat::Online<Real, Real> > stat_sigma2_;
    /** Array of the dim statistics */
    Array1D< Stat::Online<Real, Real> > stat_dim_;

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

    /** @return the standard deviation of the kth cluster */
    inline Real const& sigma2(int k) const { return sigma2_[k];}
    /** @return the dimension of the kth cluster */
    inline Real const& dim(int k) const { return dim_[k];}

    /** update statistics of the parameters. */
    void updateStatistics();
    /** Set the computed statistics */
    void setStatistics();
    /** Release the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the sigma2 and
     *  dim parameters on consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(ExprBase<Array> const& params)
    {
      for(int k=dim_.begin(); k<dim_.end(); ++k)
      { sigma2_[k] = params(k, baseIdx)   ;
        dim_[k]    = params(k, baseIdx+1);
      }
    }
};

} // namespace STK

#endif /* STK_KERNELPARAMETERS_H */
