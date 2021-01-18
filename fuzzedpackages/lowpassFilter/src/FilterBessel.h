#ifndef LOWPASSFILTER_H_FILTERBESSEL
#define LOWPASSFILTER_H_FILTERBESSEL

#include "Filter.h"
#include <Rcpp.h>
using namespace Rcpp;

/***************
* class FilterBessel
* implements abstract class Filter for the Bessel filter
* Florian Pein, 2016
***************/
class FilterBessel : public Filter {
  private:
    double truncation, C, timescaling, A;
    NumericVector a, b, c, d;
    
  public:  
    FilterBessel(const List &input);
    
    double antiderivative(const double &t) const;
 };

#endif
