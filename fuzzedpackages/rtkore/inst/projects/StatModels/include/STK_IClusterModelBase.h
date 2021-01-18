/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff

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
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_IClusterModelBase.h
 *  @brief In this file we define the abstract base class for clustering
 *  models.
 **/

#ifndef STK_ICLUSTERMODELBASE_H
#define STK_ICLUSTERMODELBASE_H

#include "STK_IStatModelBase.h"

#include <Arrays/include/STK_CArrayPoint.h>
#include <Arrays/include/STK_CArrayVector.h>

namespace STK
{
/** @ingroup StatModels
 *  @brief Base class for clustering models.
 *
 * Cluster analysis or clustering is the task of grouping a set of objects in
 * such a way that objects in the same group (called a cluster) are more similar
 * (in some sense or another) to each other than to those in other groups
 * (clusters). It is a main task of exploratory data mining, and a common
 * technique for statistical data analysis, used in many fields, including
 * machine learning, pattern recognition, image analysis, information retrieval,
 * bioinformatics, data compression, and computer graphics.
 *
 * Cluster analysis itself is not one specific algorithm, but the general task
 * to be solved. It can be achieved by various algorithms that differ
 * significantly in their notion of what constitutes a cluster and how to
 * efficiently find them. Popular notions of clusters include groups with small
 * distances among the cluster members, dense areas of the data space, intervals
 * or particular statistical distributions. Clustering can therefore be
 * formulated as a multi-objective optimization problem. The appropriate
 * clustering algorithm and parameter settings (including values such as the
 * distance function to use, a density threshold or the number of expected
 * clusters) depend on the individual data set and intended use of the results.
 * Cluster analysis as such is not an automatic task, but an iterative process
 * of knowledge discovery or interactive multi-objective optimization that
 * involves trial and failure. It is often necessary to modify data
 * preprocessing and model parameters until the result achieves the desired
 * properties.
 *
 */
class IClusterModelBase: public IStatModelBase
{
  protected:
    /** Constructor.
     * @param nbCluster number of clusters
     **/
    IClusterModelBase( int nbCluster);
    /** Constructor.
     * @param nbCluster,nbSample number of clusters and samples
     **/
    IClusterModelBase( int nbSample, int nbCluster);
    /** copy constructor.
     *  @param model the model to clone
     **/
    IClusterModelBase( IClusterModelBase const& model);

  public:
    /** destructor */
    virtual ~IClusterModelBase();

    /** @return the number of cluster */
    inline int nbCluster() const { return nbCluster_;}
    /** @return the proportions of each mixtures */
    inline CPointX const& pk() const { return pk_;};
    /** @return the zi class label */
    inline CVectorXi const& zi() const { return zi_;};

  protected:
    /** set the number of cluster of the model
     *  @param nbCluster number of cluster of the model
     * */
    inline void setNbCluster( int nbCluster) { nbCluster_ = nbCluster;}
    /** number of cluster. */
    int nbCluster_;
    /** The proportions of each mixtures */
    CPointX pk_;
    /** The zi class label */
    CVectorXi zi_;
};



} // namespace STK

#endif /* STK_ICLUSTERMODELBASE_H */


