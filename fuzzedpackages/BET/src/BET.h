/*
 * BET.h
 *
 *  Created on: Oct 20, 2020
 *      Author: wanzh
 */

#ifndef BET_H_
#define BET_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <Rcpp.h>

namespace N
{
	class BETfunction
	{
	   public:
		int d;
		int getStats();
		double getPvalue();
		std::vector<std::string> getInteraction();
		std::vector<std::vector<std::string>> getSymmInteraction();
		std::vector<int> getSymmStats();
		std::vector<std::string> getBinary();
		BETfunction(Rcpp::NumericMatrix& X_R, int depth, bool unif, bool asymptotic);

	   private:
		bool unifMargin;
    bool asympt;
		size_t numThread = 4;
		size_t p;
		std::vector<std::vector<double>> X;
		std::vector<std::string> binary_inter;
		// std::vector<std::string> symminter;
		std::vector<std::vector<std::string>> out_symminter;
		std::vector<int> out_symmstats;
		std::vector<long long> bid;
		int Stats = 0;
		double pvalue = 0;
		std::vector<std::vector<double>> imp(Rcpp::NumericMatrix& X);
		std::vector<std::string> interaction_str;
		std::vector<std::vector<int>> interactions();
		std::vector<std::string> interaction_index(bool binary);
		std::string binaryToNum(std::string binary);
		std::vector<long long> BIDs();
		std::vector<std::vector<int>> ecdf_loc(std::vector<std::vector<double>>& X);
		std::map<std::vector<int>, int> groupC(std::vector<std::vector<int>>& c);
		int locate(int c, long long b);
		std::vector<std::vector<std::vector<int>>> CBIDs(std::map<std::vector<int>, int>& count);
		std::vector<std::vector<size_t>> allComb();
		double logHypergeometricProb(double* logFacs, int a, int b, int c, int d);
		double fisherExact(int a, int b, int c, int d);
		double binomial(int n, int k, double p);
		double pbinom(int n, int k, double p);
	};
}



#endif /* BET_H_ */
