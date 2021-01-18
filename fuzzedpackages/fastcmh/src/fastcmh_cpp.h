#ifndef FASTCMH_CPP_H
#define FASTCMH_CPP_H


/* LIBRARY INCLUDES */
#include<math.h>
#include<stdlib.h>
#include<stdio.h>

//for converting int to string WITH C++11
#include <string>

// #include<Rdefines.h>
// #include<String>
#include <exception>
#include <stdexcept>

//from time_keeping.c
#include <time.h> //Already included in original LCM source code above

//cpp time
#include <ctime>

#include <R.h>
#include "Rmath.h"
#include<Rcpp.h>
#include<Rinternals.h>

//filter intervals
#include "filterIntervals.h"
#include "rcppdatawrap.h"

#include "fdr.h"

//NEED VECTOR
#include<vector>




/* CODE DEPENDENCIES */
//#include"time_keeping.c"
//#include "./chi2.h"



//this is the main function, and is exported in fastcmh_cpp.cpp
Rcpp::List main_fastcmh2(Rcpp::String, 
                         Rcpp::String, 
                         Rcpp::String, 
                         Rcpp::NumericVector,
                         Rcpp::NumericVector,
                         Rcpp::LogicalVector,
                         Rcpp::LogicalVector,
                         Rcpp::LogicalVector, 
                         Rcpp::LogicalVector);


#endif
