#include "FilterBessel.h"

/***************
 * class FilterBessel
 * implements abstract class Filter for the Bessel filter
 * Florian Pein, 2016
 ***************/
FilterBessel::FilterBessel(const List &input) : truncation(input["truncation"]), C(input["C"]),
  timescaling(input["timescaling"]), A(input["A"]), a(wrap(input["a"])), b(wrap(input["b"])),
  c(wrap(input["c"])), d(wrap(input["d"])) {}

double FilterBessel::antiderivative(const double &t) const {
  if (t <= 0.0) {
    return 0.0;
  }
  
  if (t >= truncation) {
    return 1.0;
  }
  
  double s = t * timescaling;
  NumericVector cosbs = cos(b * s);
  NumericVector sinbs = sin(b * s);
  
  return C * sum(exp(a * s) / (a * a + b * b) * (c * (a * cosbs + b * sinbs) - d * (a * sinbs - b * cosbs))) + A;
}
