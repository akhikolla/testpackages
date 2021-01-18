//=================================
// include guard
#pragma once

//=================================
// included dependencies
#include "vHMM.h"
#include "HMM.h"
#include "HMMpoisson.h"
#include "MultiGHMM.h"

/*
* note : RcppExport is an alias to ‘extern "C"‘ defined by Rcpp.
*
* It gives C calling convention to all the functions so that
* they can be called from .Call in R. Otherwise, the C++ compiler mangles the
* name of the function and .Call can’t find it.
*
* It is only useful to use RcppExport when the function is intended to be called
* by .Call. See http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
* on Rcpp-devel for a misuse of RcppExport
*/

//=================================
// wrapper functions
//=================================

//  hidden Markov model
//  Model initialiation functions
RcppExport SEXP initHMM(SEXP n, SEXP m);        //  Categorical observations
RcppExport SEXP initGHMM(SEXP n, SEXP m);       //  Univariate m = 1 or Multivariate m > 1 Continuous observations
RcppExport SEXP initPHMM(SEXP n);               //  Discrete observations
//  If the user does the model object, it verifies if it fulfills all the requirements
RcppExport SEXP verifyModel(SEXP model);
//  After the model is established, if it is necessary to change some parameters the user can use these functions
RcppExport SEXP setNames(SEXP hmm, SEXP names);
RcppExport SEXP setParameters(SEXP hmm, SEXP params);
//  Sequence evaluation
RcppExport SEXP evaluation(SEXP hmm, SEXP sequence, SEXP method);
//  Sequence decodification
RcppExport SEXP viterbi(SEXP hmm, SEXP sequence);
RcppExport SEXP forwardBackward(SEXP hmm, SEXP sequence);
//  ALearning
RcppExport SEXP loglikelihood(SEXP hmm, SEXP sequences);
RcppExport SEXP learnEM(SEXP hmm, SEXP sequences, SEXP iter, SEXP delta, SEXP pseudo, SEXP print);
RcppExport SEXP learnBW(SEXP hmm, SEXP sequences, SEXP iter, SEXP delta, SEXP pseudo, SEXP print);
//  Simulation
RcppExport SEXP generateObservations(SEXP hmm, SEXP length);
