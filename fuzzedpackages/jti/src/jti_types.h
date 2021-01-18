#ifndef JTI_TYPES_H
#define JTI_TYPES_H

#include <numeric>       // For: transform and accumulate
#include <vector>        // For: 
#include <string>        // For:
#include <algorithm>     // For: sort, std set operations etc. 
#include <math.h>        // For: log
#include <map>           // For: count_unique
#include <unordered_map> // For: dfs, mcs, count_unique
#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

using edge_pair = std::pair<int, std::pair<int, int>>;
using VI  = std::vector<int>;
using VS  = std::vector<std::string>;
using VVS = std::vector<std::vector<std::string>>;
using RL  = Rcpp::List;
using RIV = Rcpp::IntegerVector;
using RCV = Rcpp::CharacterVector;
using RCM = Rcpp::CharacterMatrix;
using RE  = Rcpp::Environment;

#endif
