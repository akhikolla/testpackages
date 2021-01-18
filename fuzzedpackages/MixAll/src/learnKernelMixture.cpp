/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

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
 * created on: 28 juil. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file learnMixture.cpp
 *  @brief In this file we launch the computation for learning a mixture model.
 **/


#include "RTKpp.h"
#include "MixAll.h"
#include <MixAll/LearnLauncher.h>

/* @param model a ClusterDiagModel S4 class
 *  @param nbCluster a vector with the number of learns to test
 *  @param models a vector of string with the model names to try
 */
RcppExport SEXP learnKmm( SEXP model, SEXP models, SEXP algo, SEXP nbCore )
{
  BEGIN_RCPP

#ifdef _OPENMP
  int cores = Rcpp::as<int>(nbCore);
  if (cores > 1) { omp_set_num_threads(cores);}
  else { omp_set_num_threads(1);}
#endif

  Rcpp::S4 s4_model(model);
  Rcpp::S4 s4_component = s4_model.slot("component");
  // build Gram matrix
  Rcpp::CharacterVector r_kernelName = s4_model.slot("kernelName");
  Rcpp::DoubleVector    r_kernelParameters = s4_model.slot("kernelParameters");
  if (!Rcpp::as<bool>(computeGramMatrix(s4_component, r_kernelName, r_kernelParameters)))
  { return Rcpp::wrap(false);}
  // create a launcher
  //STK::LearnLauncher launcher(model, models, algo);
  // return result
  //return Rcpp::wrap(launcher.run());
  return Rcpp::wrap(false);

  END_RCPP
}

