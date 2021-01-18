/**
 * @authors Hmamouche Youssef
   This file is part of NlinTS. NlinTS is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.
   NlinTS is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 **/

#ifndef NLINCAUSALITY_H
#define NLINCAUSALITY_H

#include <Rcpp.h>
#include "struct.h"
#include "varnn.h"


/*****************  Non linear causality test  *********************/
class nlinCausalityTest {
  private:

      Struct::CVDouble ts1;
      Struct::CVDouble ts2;
      double Ftest;
      unsigned lag;
      double p_value;
      double GCI;
      double criticTest;
      bool bias;

      std::vector<unsigned long> sizeOfLayersModel1;
      std::vector<unsigned long> sizeOfLayersModel2;
      VARNN univariateModel;
      VARNN bivariateModel;

  public:

      /**
      @param ts1_ the first univariate time series as a vector.
      @param ts2_ the second time series.
      @param lag_ the lag parameter.
      @param hiddenLayersOfUnivModel vector of hidden layers sizes for the univariate model.
      @param hiddenLayersOfUnivModel vector of hidden layers sizes for the bivariate model.
      @param iterations the number of iterations.
      @param learningRateInit learning rate for sgf algoritm, or initial  rate if Adam algorithm is used.
      @param bias a boolean value for the possibility of using the bias in the network.
      */

      nlinCausalityTest  (unsigned lag_ = 1);
      ~nlinCausalityTest (){};

      void buildModels (Rcpp::IntegerVector  hiddenLayersOfUnivModel,
                        Rcpp::IntegerVector  hiddenLayersOfBivModel,
                        Rcpp::StringVector activationsUnivModel = {},
                        Rcpp::StringVector activationsBivModel = {},
                        double learningRateInit = 0.1,
                        string algo = "sgd",
                        bool bias = 1);


      void fit (Rcpp::NumericVector  ts1_, Rcpp::NumericVector  ts2_, unsigned iterations);
      // The causality index
      double get_gci ();

      // Get the p-value of the test
      double get_p_value ();

      // Get the  statistic of the test
      double get_F_test () ;

      // The Summary function
      void summary ();
      };
#endif
