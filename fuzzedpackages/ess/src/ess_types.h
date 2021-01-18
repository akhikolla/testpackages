#include <numeric>       // For: 
#include <vector>        // For: 
#include <string>        // For:
#include <algorithm>     // For: sort, std set operations etc. 
#include <stack>         // For: dfs procedure
#include <math.h>        // For: log
#include <unordered_map> // For: dfs 
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

using VS  = std::vector<std::string>;
using VD  = std::vector<double>;
using VVS = std::vector<std::vector<std::string>>;
using RL  = Rcpp::List;
using RIV = Rcpp::IntegerVector;
using RCV = Rcpp::CharacterVector;
using RCM = Rcpp::CharacterMatrix;
