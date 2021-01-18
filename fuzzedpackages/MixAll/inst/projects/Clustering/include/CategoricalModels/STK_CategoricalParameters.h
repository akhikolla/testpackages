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

/** @file  STK_CategoricalParameters.h
 *  @brief In this file we define the Parameters classes for categorical
 *  mixture models
 **/

#ifndef STK_CATEGORICALPARAMETERS_H
#define STK_CATEGORICALPARAMETERS_H

#include "../STK_Clust_Util.h"

#include <Arrays/include/STK_Array1D.h>
#include <Arrays/include/STK_CArray.h>
#include <Arrays/include/STK_CArrayVector.h>

#include <STatistiK/include/STK_Stat_Online.h>

namespace STK
{

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Categorical_pjk model.
 */
template<>
struct ModelParameters<Clust::Categorical_pjk_>
{
    /** array of size nbCluster with the probabilities of the variables */
    Array1D<CArrayXX> proba_;
    /** Array of size nbCluster with the statistics of the probabilities */
    Array1D<  Stat::Online<CArrayXX, Real>  > stat_proba_;
    /** default constructor
     *  @param nbCluster number of class of the mixture
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

    // getters
    /** @return the probability of the kth cluster, jth variable, lth modality */
    inline Real const& proba(int k, int j, int l) const { return proba_[k](l,j);}
    /** @return the probabilities of the kth cluster for the jth variable */
    inline CVectorX proba(int k, int j) const { return proba_[k].col(j);}

    /** resize the set of parameter
     *  @param rangeModalities range of the  modalities
     *  @param rangeCols range of the variables
     **/
    void resize(Range const& rangeModalities, Range const& rangeCols);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** set and release the computed statistics */
    void setStatistics();
    /** Set the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the proabiliteis
     *  on consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters( ExprBase<Array> const& params)
    {
      int kp = params.beginRows();
      for(int k=proba_.begin(); k<proba_.end(); ++k)
      {
        for (int l = proba_[k].beginRows(); l < proba_[k].endRows(); ++l, ++kp)
        {
          for (int j = proba_[k].beginCols(); j < proba_[k].endCols(); ++j)
          { proba_[k](l, j) = params(kp , j) ;}
        }
      }
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a Categorical_pk model.
 */
template<>
struct ModelParameters<Clust::Categorical_pk_>
{
    /** array of size nbCluster with the probabilities of the variables */
    Array1D<CVectorX> proba_;
    /** array of size nbCluster with the statistics of the probabilities */
    Array1D<  Stat::Online<CVectorX, Real>  > stat_proba_;

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

    // getters
    /** @return the probability of the kth cluster, jth variable, lth modality */
    inline Real const& proba(int k, int j, int l) const { return proba_[k][l];}
    /** @return the probabilities of the kth cluster for the jth variable */
    inline CVectorX proba(int k, int j) const { return proba_[k];}

    /** resize the set of parameter
     *  @param rangeModalities range of the  modalities
     *  @param rangeCols range of the variables (not used)
     **/
    void resize(Range const& rangeModalities, Range const& rangeCols);

    /** update statistics of the parameters. */
    void updateStatistics();
    /** set and release the computed statistics */
    void setStatistics();
    /** Set the computed statistics */
    void releaseStatistics();

    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters( ExprBase<Array> const& params)
    {
      int kp = params.beginRows();
      for(int k=proba_.begin(); k<proba_.end(); ++k)
      {
        for (int l = proba_[k].beginRows(); l < proba_[k].endRows(); ++l, ++kp)
        {
          proba_[k][l] = 0.;
          for (int j = params.beginCols(); j < params.endCols(); ++j)
          { proba_[k][l] += params(kp , j) ;}
          proba_[k][l] /= proba_[k].sizeCols();
        }
      }
    }
};


} // namespace STK

#endif /* STK_CATEGORICALPARAMETERS_H */
