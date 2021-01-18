#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false
#include <Rcpp.h>
#include "wdm.hpp"

// [[Rcpp::export]]
double wdm_cpp(const std::vector<double>& x,
               const std::vector<double>& y,
               std::string method,
               const std::vector<double>& weights,
               bool remove_missing)
{
    return wdm::wdm(x, y, method, weights, remove_missing);
}

std::vector<double> convert_vec(const Rcpp::NumericVector& x)
{
    return Rcpp::as<std::vector<double>>(x);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix wdm_mat_cpp(const Rcpp::NumericMatrix& x,
                                std::string method,
                                const std::vector<double>& weights,
                                bool remove_missing)
{
    using namespace Rcpp;
    size_t d = x.ncol();
    if (d == 1)
        throw std::runtime_error("x must have at least 2 columns.");

    NumericMatrix ms(d, d);
    for (size_t i = 0; i < x.cols(); i++) {
        for (size_t j = i; j < x.cols(); j++) {
            if (j == i) {
                ms(j, i) = 1.0;
                continue;
            }
            ms(i, j) = wdm::wdm(convert_vec(x(_, i)),
                                convert_vec(x(_, j)),
                                method,
                                weights,
                                remove_missing);
            ms(j, i) = ms(i, j);
        }
    }

    return ms;
}

// [[Rcpp::export]]
Rcpp::List indep_test_cpp(const std::vector<double>& x,
                          const std::vector<double>& y,
                          std::string method,
                          const std::vector<double>& weights,
                          bool remove_missing,
                          std::string alternative)
{
    wdm::Indep_test test(x, y, method, weights, remove_missing, alternative);
    return Rcpp::List::create(
        Rcpp::Named("estimate") = test.estimate(),
        Rcpp::Named("statistic") = test.statistic(),
        Rcpp::Named("p_value") = test.p_value(),
        Rcpp::Named("n_eff") = test.n_eff(),
        Rcpp::Named("method") = method,
        Rcpp::Named("alternative") = alternative
    );
}
