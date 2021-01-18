#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
using namespace std;

//' Calculates DysPIA statistic values for the gene pair sets
//'
//' @name calcDyspiaStatCumulativeBatch
//' @param stats Named numeric vector with gene pair-level statistics 
//'        sorted in decreasing order (order is not checked).
//' @param DyspiaParam DysPIA weight parameter (0 is unweighted, suggested value is 1).
//' @param pathwayScores Vector with enrichment scores for the pathways in the database.
//' @param pathwaysSizes Vector of pathway sizes.
//' @param iterations Number of iterations.
//' @param seed Seed vector
//' @return List of DysPIA statistics for gene pair sets.
//' @export
// [[Rcpp::export]]
List calcDyspiaStatCumulativeBatch(
        NumericVector const& stats,
        double DyspiaParam,
        NumericVector const& pathwayScores,
        IntegerVector const& pathwaysSizes,
        int iterations,
        int seed);


//' Calculates DysPIA statistic values for all the prefixes of a gene pair set
//'
//' @name calcDyspiaStatCumulative
//' @param stats Named numeric vector with gene pair-level statistics
//'        sorted in decreasing order (order is not checked)
//' @param selectedStats indexes of selected gene pairs in a `stats` array
//' @param DyspiaParam DysPIA weight parameter (0 is unweighted, suggested value is 1)
//' @return Numeric vector of DysPIA statistics for all prefixes of selectedStats.
//' @export
// [[Rcpp::export]]
NumericVector calcDyspiaStatCumulative(
        NumericVector const& stats,
        IntegerVector const& selectedStats, // Indexes start from one!
        double DyspiaParam
);

