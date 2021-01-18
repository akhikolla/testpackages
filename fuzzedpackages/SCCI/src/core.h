#ifndef CORE_H
#define CORE_H

#include <Rcpp.h>
#include <vector>
#include <set>
#include <assert.h>

/**
 * Computes log2 with log2(0) = 0 (used in entropy calculation for 0 * log2(0))
 **/
double myLog2(double);

/**
 * Helpers
 **/
double log2fac(int);
double log2nChoosek(int,int);
double stirling(int,int);

/**
 * Helper functions to compute the costs for the stochastic complexity
 **/
double binaryRegretPrecal(int);
double regretPrecal(int,int);

/**
 * Calculates the regred term of the stochastic complexity with
 * @M: number of elements
 * @K: number of categories
 **/
double regret(int,int);

/**
 * Computes the entropy given:
 * @counts: a map between value and its count
 * @sum: total number of counts
 **/
double entropy(std::map<int,int>&,int);
double conditionalEntropy(const std::vector<int>&,const std::vector<int>&);

/**
 * Computes the entropy given:
 * @counts: a vector of counts followed by the sum of all counts in the last cell
 * Important note: the last cell contains the sum of all counts
 **/
double entropy(std::vector<int>&);

/**
 * Computes the stoachastic complexity of a vector.
 **/
double SC(Rcpp::IntegerVector&);
double SC(std::vector<int>&);

/**
 * Computes the conditional stoachastic complexity between two vectors.
 **/
double conditionalSC(std::vector<int>&,std::vector<int>&);
double conditionalSC(Rcpp::IntegerVector&,std::vector<int>&);

/**
 * Input: vector of categories
 * Output: vector of categories starting from 0 and going to domain - 1;
 *         THE LAST CELL CONTAINS THE DOMAIN SIZE
 * Important Note: This method does not consider an order in the data,
 *                 but is only for unordered categorical data--no ordinal
 **/
std::vector<int> getNiceCategories(std::vector<int>&);
std::vector<int> getNiceCategories(Rcpp::IntegerVector&);

#endif // CORE_H
