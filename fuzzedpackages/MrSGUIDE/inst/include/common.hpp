//
// Created by Peigen Zhou on 8/1/18.
//

#ifndef SUBGUIDE_SUBGUIDE_H
#define SUBGUIDE_SUBGUIDE_H

#define ARMA_USE_CXX11

#include <iostream>
#include "RcppArmadillo.h"
#include <cmath>
#include <vector>
#include <string>

const int miss = 12345679;
// const int seed = 1234;

// arma::arma_rng::set_seed(seed);

#define cout Rcpp::Rcout
#define cerr Rcpp::Rcerr

using std::endl;

#ifdef _OPENMP
#include <omp.h>
#endif

#define BOOST_DISABLE_ASSERTS
/*
#include "../inst/include/spdlog/spdlog.h"
#include "../inst/include/spdlog/fmt/ostr.h"
#include "../inst/include/spdlog/sinks/basic_file_sink.h"
#include "../inst/include/spdlog/sinks/stdout_sinks.h"
*/

#endif //SUBGUIDE_SUBGUIDE_H
