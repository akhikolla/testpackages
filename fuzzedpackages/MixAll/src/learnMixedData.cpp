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

/** @file learnMixedData.cpp
 *  @brief In this file we launch the computation for learning a mixture model with mixed data.
 **/


#include "RTKpp.h"
#include "MixAll.h"
#include <MixAll/LearnLauncher.h>


/*
 *  @param model ClusterMixedDataModel S4 class
 *  @param critName name criteria string
 *  @param nbCore number of core to use
 */
RcppExport SEXP learnMixedData( SEXP model, SEXP algo, SEXP nbCore  )
{
  BEGIN_RCPP

#ifdef _OPENMP
  int cores = Rcpp::as<int>(nbCore);
  if (cores > 1) { omp_set_num_threads(cores);}
  else { omp_set_num_threads(1);}
#endif
  Rcpp::S4 s4_model(model);
  Rcpp::S4 s4_algo(algo);
  // create a launcher
  STK::LearnLauncher launcher(s4_model, s4_algo);
  // return result
  return Rcpp::wrap(launcher.run());

  END_RCPP
}
