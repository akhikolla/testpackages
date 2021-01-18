/*!
 *
 * @author youssef hmamouche
 *
   This file is part of NlinTS. NlinTS is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.
   NlinTS is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   *
 */

#pragma once


#include <Rcpp.h>
#include "struct.h"
#include "varnn.h"

class VARNN_Export
{
private:
    Struct::CMatDouble M;
    VARNN Obj;

public:
    VARNN_Export (
    			        unsigned p,
				          Rcpp::IntegerVector,
                  Rcpp::StringVector activations,
                  double learning_rate_init,
                  string  algo,
                  bool bias);

    VARNN_Export(){};
   ~VARNN_Export(){};
    Rcpp::DataFrame forecast (Rcpp::DataFrame DF);
    void fit (Rcpp::DataFrame, int);
    void train (Rcpp::DataFrame DF);
    Rcpp::NumericVector getSSR ();
};
