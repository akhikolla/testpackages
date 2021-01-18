#ifndef RCPPDATAWRAP_H
#define RCPPDATAWRAP_H

#include<vector>
#include<numeric>
#include <R.h>
#include<Rcpp.h>
#include<Rinternals.h>
#include "filterIntervals.h"

using std::vector;
using std::size_t;

//extract the individual vectors into a dataframe
Rcpp::DataFrame extractDataFrameFromIntervalVector(const vector<Interval>& );

//create an empty data frame, in case we need to return it
Rcpp::DataFrame createEmptyDataFrame();


//just to wrap tau, l, pvalue
Rcpp::DataFrame createDataFrameTauLPvalue(const vector<long long>& , const vector<long long>& , const vector<double>& );


//create a return list that indicates an error occurred
Rcpp::DataFrame createErrorReturnList();


//used for testing the filtering
//this function is exported to R
Rcpp::DataFrame cpp_test_filtering(const Rcpp::DataFrame&);


#endif
