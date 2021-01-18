// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppProgress)]]

#include <progress.hpp>
#include <Rcpp.h>
#include "clock.h"

using namespace Rcpp;

// [[Rcpp::export]]
long double do_timing(const RObject& expr, const Environment& env) {
    time_point start = Clock::now();
    Rf_eval(expr, env);
    time_point end = Clock::now();
    return duration_cast<Seconds>(end - start).count();
}

// [[Rcpp::export]]
NumericVector do_benchmark(const List& exprs, const Environment& env, const IntegerVector& order,
                           bool gc = false, bool progress = false) {
    std::size_t n = order.length();
    NumericVector res = no_init(n);
    Progress pb(n, progress);
    for (std::size_t i = 0; i < n; ++i) {
        pb.increment();
        RObject expr = exprs[order[i] - 1];
        res[i] = do_timing(expr, env);
        Progress::check_abort();
        if (gc) R_gc();
    }
    return res;
}
