// [[Rcpp::plugins(cpp11)]]

#include "clock.h"

//' Get timer precision.
//' @return Number of seconds from one clock tick to the next.
//' @export
// [[Rcpp::export]]
double timer_precision() {
    return static_cast<double>(Clock::period::num) / Clock::period::den;
}
