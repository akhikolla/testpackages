/**
 * @authors Hmamouche Youssef
 *
 This file is part of NlinTS. NlinTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 2 of the License, or
 (at your option) any later version.
 NlinTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 **/

#include <Rcpp.h>
#include "../inst/include/nlinCausalityTest.h"
#include "../inst/include/Ftable.h"
#include "../inst/include/fdist.h"
#include "../inst/include/varnn.h"

using namespace Struct;
using namespace std;

/************************************/
nlinCausalityTest::nlinCausalityTest  (unsigned lag_) {
  // Checking if the lag is positif
    if (lag_ <= 0)
        throw string ("Error: The lag value is incorrect, try a strictly positive value.");

    lag  = lag_;
}

/************************************/
void nlinCausalityTest::buildModels ( Rcpp::IntegerVector  hiddenLayersOfUnivModel,
                                      Rcpp::IntegerVector  hiddenLayersOfBivModel,
                                      Rcpp::StringVector activationsUnivModel,
                                      Rcpp::StringVector activationsBivModel,
                                      double learningRateInit,
                                      string algo,
                                      bool bias_)
{

    bias = bias_;
    // Size of hidden layers
    sizeOfLayersModel1 = Rcpp::as < vector<unsigned long>  > (hiddenLayersOfUnivModel);
    sizeOfLayersModel2 = Rcpp::as < vector<unsigned long>  > (hiddenLayersOfBivModel);

    // Activations vectors
    vector<string> activationsModel1 = Rcpp::as <vector<string> > (activationsUnivModel);
	  vector<string> activationsModel2 = Rcpp::as <vector<string> > (activationsBivModel);

    univariateModel = VARNN (sizeOfLayersModel1, lag, bias, learningRateInit, activationsModel1, algo);
    bivariateModel = VARNN (sizeOfLayersModel2, lag, bias, learningRateInit, activationsModel2, algo);
}

/************************************/
void nlinCausalityTest::fit (Rcpp::NumericVector  ts1_, Rcpp::NumericVector  ts2_, unsigned iterations)
{
    for (const auto val:ts1_)
        ts1.push_back (val);

    for (const auto val:ts2_)
        ts2.push_back (val);

    if (ts1.size() != ts2.size())
            throw string ("Error: The variables have not the same length.");

    // variables
    unsigned nl (ts1.size ());
    double RSS1(0), RSS0(0);

    CMatDouble M;
    //M.reserve (2);
    M.push_back (ts1);

    // Fit the model on univariate ts
    univariateModel.fit (M, iterations);

    // add second variable
    M.push_back(ts2);
    // Fit Model2
    bivariateModel.fit (M,  iterations);

    // Evaluate the two models
    RSS0 = univariateModel.getSSR () [0];
    RSS1 = bivariateModel.getSSR ()  [0];

  // Granger causality index
    GCI = log (RSS0 / RSS1);

    // RSS1 > RSS0, adding the predictive variable degrade the predictions, in this case we put GCI = 0
    // In this case, we dont need a statistical test
    if (GCI < 0)
        GCI = 0;

    int T = nl - lag;

    // Compute the number of parameters of the univariate and bivariate models
    int nbParamsU = 0, nbParamsB = 0, lastLayer = lag + 1;

    for (auto val:sizeOfLayersModel1)
    {
        nbParamsU += val * lastLayer + (int) (bias);
        lastLayer = val;
    }
    nbParamsU += lastLayer;
    lastLayer = 2*lag + 1;

    for (auto val:sizeOfLayersModel2)
    {
        nbParamsB += val * lastLayer + (int) (bias);
        lastLayer = val;
    }

    nbParamsB += lastLayer;

    if (nbParamsB > T)
        nbParamsB = T;

    // degrees of liberty for the F-test
    int d1, d2;
    d1 = nbParamsB - nbParamsU;
    d2 = T - nbParamsB;

    if (d2 <= 0)
    {
        Ftest = nan ("");
        p_value = nan ("");
        criticTest = nan ("");
    }

    else
    {
        /* F test */
        Ftest = ((RSS0 - RSS1) / RSS1) *  ( (double) d2 / (double) d1 );

        /* p-value of the F-test */
         p_value = getPvalue (Ftest , d1 , d2);

        if (d1 <= 20 and d2 <= 100)
            criticTest = ftable[d2][d1];
        else if (d1 > 20 and d2 <= 100)
            criticTest = ftable[d2][20];
        else if (d1 <= 20 and d2 > 100)
            criticTest = ftable[100][d1];
        else if (d1 > 20 and d2 > 100)
            criticTest = ftable[100][20];
    }
}
void nlinCausalityTest::summary ()
{
    Rcpp::Rcout <<  "--------------------\n";
    Rcpp::Rcout <<  "    Test summary"  << "\n";
    Rcpp::Rcout <<  "--------------------\n";
    Rcpp::Rcout <<  "The lag parameter: p = "<< lag << "\n";
    Rcpp::Rcout <<  "The Granger causality Index: GCI = "<< GCI << "\n";
    Rcpp::Rcout <<  "The value of the F-test: "<< Ftest << "\n";
    Rcpp::Rcout <<  "The p_value of the F-test: "<< p_value << "\n";
    Rcpp::Rcout <<  "The critical value at 5% of risk: "<< criticTest <<"\n";
}
double nlinCausalityTest::get_p_value () {
    return p_value;
}

double nlinCausalityTest::get_gci () {
    return GCI;
}

double nlinCausalityTest::get_F_test () {
    return Ftest;
}
