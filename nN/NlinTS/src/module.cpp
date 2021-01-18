/*
  NlinTS -- Rcpp package for non-linear time series analysis
  Copyright (C) 2017 - 2020  Hmamouche youssef
  This file is part of NlinTS
  NlinTS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.
  NlinTS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include<Rcpp.h>
#include "../inst/include/varnnExport.h"
#include "../inst/include/tests.h"
#include "../inst/include/causalityTest.h"
#include"../inst/include/nlinCausalityTest.h"

using namespace Rcpp;

// 'The VARMLP module
RCPP_MODULE (VAR_MLP) {
  class_<VARNN_Export> ("VARNN_Export")
        .constructor <unsigned, Rcpp::IntegerVector, Rcpp::StringVector, double, string, bool> ()
        .property("MSE", &VARNN_Export::getSSR, "mean squared error of the model on training data")
        .method ("fit", &VARNN_Export::fit, "fit the model")
        .method ("train", &VARNN_Export::train, "Update the model from input data")
        .method ("forecast", &VARNN_Export::forecast, "Computes the predictions");
}

// 'The dickey fuller test module.
RCPP_MODULE (DickeyFuller) {
  class_<DickeyFuller> ("DickeyFuller")
        .constructor <Rcpp::NumericVector, int> ()
        .method ("summary", &DickeyFuller::summary, "Summary of the test")
        .property ("df", &DickeyFuller::getDF, "return the value of test");
}

// 'The Granger causality test module.
RCPP_MODULE (causalityTest) {
  class_<causalityTest> ("causalityTest")
        .constructor <Rcpp::NumericVector, Rcpp::NumericVector, int, bool> ()
        .method ("summary", &causalityTest::summary, "Summary of the test")
        .property ("pvalue", &causalityTest::get_p_value, "return the p-value of the test")
        .property ("gci", &causalityTest::get_gci, "return the granger causality index of the test")
        .property ("Ftest", &causalityTest::get_F_test, "return the value of F test");
}


//'The non linear Granger CausalityTest module.
RCPP_MODULE (nlinCausalityTest) {
    class_<nlinCausalityTest> ("nlinCausalityTest")
        .constructor <unsigned> ()
        .method ("buildModels", &nlinCausalityTest::buildModels, "Build the two models.")
        .method ("fit", &nlinCausalityTest::fit, "fit the test")
        .method ("summary", &nlinCausalityTest::summary, "Summary of the test")
        .property ("gci", &nlinCausalityTest::get_gci, "return the granger causality index of the test")
        .property ("pvalue", &nlinCausalityTest::get_p_value, "returns the p-value of the test")
        .property ("Ftest", &nlinCausalityTest::get_F_test, "returns the value of F test");
}
