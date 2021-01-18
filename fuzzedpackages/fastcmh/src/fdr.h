#ifndef FDR_TARONE_H
#define FDR_TARONE_H

//based on Gilbert (2005) A modified false discovery rate multiple-comparisons procedure for discrete data, applied to human immunodeficiency virus genetics

//Dean Bodenham June 2016

#include<iostream>
#include<vector>
//for sort
#include<algorithm>
#include<math.h>

//NOTE:
//'difficult' part was sorting the three vectors pvalue, tau and l according the the order in pvalue. Although it would have been relatively straightforward to use the sort function and a lambda function, for some reason the g++ compiler on MAc OSX (gcc 4.2) would not accept this. Then, after installing g++ 4.5, it created other errors. In the end, I decided to make a pair of pvalue and a perm vector (could have used iota), sort the pair, and then extract the permutation half of the pair. Later, access elements of the three vectors according to the permutation order.
//Based on method from: 
// http://stackoverflow.com/questions/236172/how-do-i-sort-a-stdvector-by-the-values-of-a-different-stdvector


std::vector<long long> gilbertFDR(std::vector<double>& pvalue, std::vector<long long>& tau, std::vector<long long>& l, double alpha, bool useDependence);

std::vector<double> extractFdrPvalue(const std::vector<double>& pvalue, const std::vector<long long>& perm);

std::vector<long long> extractFdrTau(const std::vector<long long>& tau, const std::vector<long long>& perm);

std::vector<long long> extractFdrL(const std::vector<long long>& l, const std::vector<long long>& perm);

#endif
