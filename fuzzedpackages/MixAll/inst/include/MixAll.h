/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/
/*
 * Project:  MixAll
 * created on: 28 juil. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/
/** @file MixAll.h
 *  @brief header file with all methods that can be used with a .Call from
 *  the package MixAll
 **/

#ifndef MIXALL_H
#define MIXALL_H

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

// mixture models
/** @param model any model derived from IClusterModel S4 class
 *  @param nbCluster a vector with the number of clusters to test
 *  @param models a vector of string with the model names to test
 *  @param nbCore number of core to use
 */
SEXP clusterMixture( SEXP model, SEXP nbCluster, SEXP models, SEXP nbCore);
/** @param model ClusterMixedDataModel S4 class
 *  @param nbCluster a vector with the number of clusters to test
 *  @param nbCore number of core to use
 */
SEXP clusterMixedData( SEXP model, SEXP nbCluster, SEXP nbCore);

/** @param model  any model derived from IClusterModel S4 class
 *  @param result an instance of ClusterPredict S4 class
 *  @param nbCore number of core to use
 */
SEXP clusterPredict( SEXP model, SEXP result, SEXP nbCore );

// kmm methods
/** estimate the kernel mixture model
 *  @param model a Kmm S4 class
 *  @param nbCluster a vector with the number of clusters to test
 *  @param models a vector of string with the model names to try
 *  @param critName selection model criteria (string)
 *  @param nbCore number of core to use
 */
SEXP kmm( SEXP model, SEXP nbCluster, SEXP models, SEXP nbCore );
/** estimate the kernel mixture model for mixed data
 *  @param model a Kmm S4 class
 *  @param nbCluster a vector with the number of clusters to test
 *  @param nbCore number of core to use
 */
SEXP kmmMixedData( SEXP model, SEXP nbCluster, SEXP nbCore  );
/** Compute the Gram matrix and overwrite data slot in S4 component with the result
 *  @param component a KmmModelComponent S4 class
 *  @param kernelName a string with the name of the kernel to use
 *  @param kernelParameters a vector with the optional parameters
 */
SEXP computeGramMatrix( SEXP component, SEXP kernelName, SEXP kernelParameters);

// learning methods
/** @param model ClusterDiagModel S4 class
 *  @param models a vector of string with the model names to try
 *  @param critName selection model criteria (string)
 *  @param nbCore number of core to use
 */
SEXP learnMixture( SEXP model, SEXP models, SEXP algo, SEXP nbCore);
/** @param model ClusterMixedDataModel S4 class
 *  @param algo estimation strategy S4 class
 *  @param critName selection model criteria (string)
 *  @param nbCore number of core to use
 */
SEXP learnMixedData( SEXP model, SEXP algo, SEXP nbCore);
/** @param model a ClusterKmm S4 class
 *  @param models a vector of string with the model names to try
 *  @param nbCore number of core to use
 */
SEXP learnKmm( SEXP model, SEXP models, SEXP algo, SEXP nbCore);
/** @param model a ClusterKmm S4 class
 *  @param models a vector of string with the model names to try
 *  @param critName name criteria string
 *  @param nbCore number of core to use
 */
//SEXP learnKmmMixedData( SEXP model, SEXP algo, SEXP nbCore);

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* MIXALL_H */
