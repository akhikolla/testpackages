#ifndef STEPR_H_DYANAMICPROGRAMS
#define STEPR_H_DYANAMICPROGRAMS

#include <Rcpp.h>

#include "Data.h"
#include "IntervalSystem.h"

using namespace Rcpp;

List fitSimpleDynamicProgram(Data * const data, IntervalSystem * const intervalSystem);
List fitIntervalDynamicProgram(Data * const data, IntervalSystem * const intervalSystem);
List fitBandDynamicProgram(Data * const data, IntervalSystem * const intervalSystem);
NumericVector computeStatistic(Data * data, const List &input);
List findSmallScales(Data * data, const List input);

#endif
