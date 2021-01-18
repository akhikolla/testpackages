/*--------------------------------------------------------------------*/
/*     Copyright (C) 2013-2014  Loic Yengo

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : loic.yengo@gmail.com
*/

/*
 * Project:  clere
 * created on: July 10, 2013
 * Author: Loic Yengo
 *
 **/

#include <RcppEigen.h>
#include "IO.h"
#include "Fit.h"
#include "ModelSelect.h"


RcppExport SEXP clere(SEXP R_IOobj){  
  BEGIN_RCPP  
    Rcpp::S4 obj(R_IOobj);
  IO io = IO(obj);
  // Set random number generator
  // int seed = io.seed;
  // std::srand(seed);
  
  if(io.instantiated!=0){
    if(io.analysis=="fit"){
      Fit Analysis = Fit(&io);
      Analysis.fitModel();
      Analysis.output();
    }
    if(io.analysis=="aic" or io.analysis=="bic" or io.analysis=="icl"){
      ModelSelect Analysis = ModelSelect(&io);
      Analysis.fitAllModels();
      Analysis.findBestModel();
      Analysis.output();
    }
    io.updateObj(obj);
  }
  //return obj;
  END_RCPP
}
