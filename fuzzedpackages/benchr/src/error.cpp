// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include "clock.h"

using namespace Rcpp;

void do_nothing() {}

// [[Rcpp::export]]
long double timer_error(std::size_t rounds = 2e5) {
    bool observed = false;
    const duration zero = duration::zero();
    const duration max = duration::max();
    duration error = max;
    for (std::size_t i = 0; i < rounds; ++i) {
        time_point start = Clock::now();
        do_nothing();
        time_point end = Clock::now();
        duration diff = end - start;
        if (diff > zero && diff < error) {
            observed = true;
            error = diff;
        }
        checkUserInterrupt();
    }
    if (!observed) {
        warning("Could not measure timer error. Your clock might lack precision.");
        error = zero;
    } else if (error == max)
        stop("Observed error too large.");
    return duration_cast<Seconds>(error).count();
}
