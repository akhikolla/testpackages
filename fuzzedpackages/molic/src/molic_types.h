#include <iostream>      // For:
#include <numeric>       // For: 
#include <vector>        // For: 
#include <string>        // For:
#include <algorithm>     // For: sort, std set operations etc. 
#include <stack>         // For: dfs procedure
#include <Rcpp.h>        // For: Interface to R
#include <regex>         // For: function na_b : the b'th slice in table a
#include <math.h>        // For: log
#include <map>           // For: count_unique
// #include <unordered_map> // For: dfs 
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

using VI  = std::vector<int>;
using VD  = std::vector<double>;
using VS  = std::vector<std::string>;
using VVS = std::vector<std::vector<std::string>>;
using MSI = std::map<std::string, int>;
using RL  = Rcpp::List;
using RNV = Rcpp::NumericVector;
using RIV = Rcpp::IntegerVector;
using RCV = Rcpp::CharacterVector;
using RCM = Rcpp::CharacterMatrix;
