//============================================================================
// Name        : BET2nd.cpp
// Author      : wanzhang
// Version     :
// Copyright   : Your copyright notice
// Description : BET 2nd version in C++, Ansi-style
//============================================================================

#include "BET.h"
#include <bitset>
#include <iostream>
#include <math.h>
#include <ctime>
#include <cstdlib>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <limits>
#include <fstream>
#include <sstream>
#include <map>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
using namespace N;

vector<vector<double> > BETfunction::imp(NumericMatrix& X)
{
  int r = X.nrow(), c = X.ncol();
  vector<vector<double>> v(r, vector<double>(c));
  for (int i = 0; i < r; i++){
    for (int j = 0; j < c; j++)
      v[i][j] = X(i, j);
  }
  return v;
}

vector<vector<int> > BETfunction::interactions()
{
	// all interactions on depth = d
	vector<vector<int>> res;
	vector<int> nums(d);
	for (int i = 0; i < d; i++){
		nums[i] = i + 1;
	}
	for(int num: nums){
	  int n = res.size();
	  res.push_back({num});

	  for (int i = 0; i < n; ++i){
		  auto temp = res[i];
		  temp.push_back(num);
		  res.push_back(temp);
	  }
	}
	return res;
}

vector<string> BETfunction::interaction_index(bool binary)
{
	// interactions in string: 1-d and binary
	vector<vector<int>> itr = interactions();
	vector<string> res;
	string s, str;
	string str0(d, '0');
	res.push_back(str0);
	for (size_t i = 0; i < itr.size(); i++)
	{
		s = accumulate(itr[i].begin()+1, itr[i].end(), to_string(itr[i][0]),
				[](const string& a, int b){return a + ':' + to_string(b);});
		// binary
		str = str0;
		for (size_t j = 0; j < itr[i].size(); j++){
			str.replace(itr[i][j]-1, 1, "1");
		}
		if (binary){
			res.push_back(str);
		}else{
			res.push_back(s + " " + str);
		}

	}
	return res;
}

vector<long long> BETfunction::BIDs()
{
	int n = (int)round(pow(2, d-1));
	vector<long long> b;
	string str1(n, '1');
	string str2(n, '0');
	//	basic configurations
	b.push_back(stoll(str2 + str1, nullptr, 2)); // A0
	for (int i = 1; i < d; i++)
	{
	//	recursively: A_i = A_{i-1} XOR (A_{i-1} << (1 << (d-i)))
		b.push_back(b[i-1] ^ (b[i-1] << (1 << (d-i-1))));
	}

	vector<vector<long long>> subs;
	//	all subsets of basic configurations
	for (long long num: b){
		int n = subs.size();
		subs.push_back({num});

		for (int i = 0; i < n; ++i){
			auto temp = subs[i];
			temp.push_back(num);
			subs.push_back(temp);
		}
	}

	double len = (int)round(pow(2, d));
	string str((int)len, '1');
	unsigned long long a0 = stoull(str, nullptr, 2);
	vector<long long> res;
	long long x;
	for (size_t i = 0; i < subs.size(); i++)
	{
	//  NOT XOR
		x = a0;
		for (size_t j = 0; j < subs[i].size(); j++){
			x = ~(x ^ subs[i][j]);
		}

		x = x & a0; // keep only 2^d digits
		res.push_back(x);
	}
	return res; // BIDs on one dimension
}

vector<vector<int>> BETfunction::ecdf_loc(vector<vector<double>>& X)
{
	int n = X.size(), p = X[0].size();
	vector<vector<int>> C(n, vector<int>(p));
	double u = 0;
	vector<double> sorted(n);
	if (p == 1)
	{
		for (int i = 0; i < n; i++){
			// throw error in R
//			if (X[i][0] > 1 or X[i][0] < 0){
//				// test whether unif[0,1]
//				throw "Data out of range [0, 1]!";
//			}
			C[i][0] = (int) ceil(X[i][0] * (int)round(pow(2, d)));
		}
	}
	else{
		if (unifMargin)
			// already uniformed
			for (int j = 0; j < p; j++){
				for (int i = 0; i < n; i++){
					C[i][j] = (int)X[i][j];
				}
			}
		else
		{
			for (int j = 0; j < p; j++){
				for (int i = 0; i < n; i++){
					sorted[i] = X[i][j];
				}
				// store sorted data
				sort(sorted.begin(), sorted.end(), greater<double>());
				for (int i = 0; i < n; i++){
					//	empirical cdf
					vector<double>::iterator ans = find(sorted.begin(), sorted.end(), X[i][j]);
					int index = distance(sorted.begin(), ans);
					u = n - index;
					C[i][j] = (int) ceil(u/n * (int)round(pow(2, d)));
				}
			}
		}

	}
	return C; // Cij
}

map<vector<int>, int> BETfunction::groupC(vector<vector<int>>& c)
{
	map<vector<int>, int> count;
	for (auto& v: c){
		count[v]++;
	}
	return count;
}

int BETfunction::locate(int c, long long b)
{
	// the c th bits of b (BID)
	int move = (int)round(pow(2, d)) - c;
	int res = (b >> move) & 1;
	return res;
}

vector<vector<vector<int>>> BETfunction::CBIDs(map<vector<int>, int>& count)
{
	size_t numCount = count.size(), numBIDs = (int)round(pow(2, d));
	// store all Cij * BIDs in to a matrix: p * #BIDs * #Cijs
	vector<vector<vector<int>>> res(p, vector<vector<int>>(numBIDs, vector<int>(numCount)));

	for (int i = 0; i < (int)p; i++){
		size_t j = 0;
		vector<vector<int>> temp(numBIDs, vector<int>(numCount, 1));
		// go over all BIDs:
		for (map<vector<int>,int>::iterator it=count.begin(); it!=count.end()&&(j < numCount); ++it, j++){
			// go over all p dimensions:
			for (size_t k = 0; k < numBIDs - 1; k++){
				temp[k+1][j] = locate(it->first[i], bid[k]);
			}
		}
		res[i] = temp;
	}
	return res;
}

vector<vector<size_t>> BETfunction::allComb()
{
	// # of BIDs
	size_t b = (int)round(pow(2, d));

	// # of combinations = 2^d^p
	size_t rows = (int)round(pow(b, p));

	// store all combinations
	vector<vector<size_t>> res(rows, vector<size_t>(p));

	vector<size_t> indices(p);
	vector<size_t> choose(p);

	// p columns, each column from 0 to 2^d-1
	vector<vector<size_t>> comb(p, vector<size_t>(b));
	for (size_t i = 0; i < p; i++){
		for (size_t j = 0; j < b; j++){
			comb[i][j] = j;
		}
	}

	size_t i = 0;
	size_t k = 0;
	while (i < p){
		for (size_t j = 0; j < p; j++){
			choose[j] = comb[j][indices[j]];
		}

		// store into res
		res[k] = choose;

		i = 0;
		while (i < comb.size() && ++indices[i] == comb[i].size())
			indices[i++] = 0;

		k++;
	}
	return res;
}

double BETfunction::logHypergeometricProb(double* logFacs, int a, int b, int c, int d)
{
	return logFacs[a+b] + logFacs[c+d] + logFacs[a+c] + logFacs[b+d] - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a+b+c+d];
}

double BETfunction::fisherExact(int a, int b, int c, int d)
{
	int n = a + b + c + d;
	double* logFacs = new double[n+1] ();

	for(int i=1; i < n+1; ++i) {
		logFacs[i] = logFacs[i-1] + log((double)i); //store factorization
	}
	double logpCutoff = logHypergeometricProb(logFacs,a,b,c,d); // *** logFacs added
	double pFraction = 0;
	for(int x=0; x <= n; ++x) {
	if ( a+b-x >= 0 && a+c-x >= 0 && d-a+x >=0 ) {
		double l = logHypergeometricProb(logFacs, x,a+b-x,a+c-x,d-a+x);
		if ( l <= logpCutoff )
			pFraction += exp(l - logpCutoff);
		}
	}
	double logpValue = logpCutoff + log(pFraction);
	double hp = exp(logpCutoff)/2;
	delete [] logFacs;
	return exp(logpValue) - hp;
}

double BETfunction::binomial(int n, int k, double p)
{
  // binomial distribution with dp
  if (n < 0 || k < 0)
    return 0;
  // store previous results
  vector<vector<double> > dp(n+1, vector<double>(k+1));
  // marginal condition
  dp[0][0] = 1;
  for (int i = 1; i < n+1; ++i){
    dp[i][0] = (1 - p) * dp[i-1][0];
  }
  for (int j = 1; j < k+1; ++j){
    dp[0][j] = 0;
  }
  for (int i = 1; i < n+1; ++i){
    for (int j = 1; j < k+1; ++j){
      // recursion
      dp[i][j] = (1 - p) * dp[i-1][j] + p * dp[i-1][j-1];
    }
  }
  return dp[n][k];
}

double BETfunction::pbinom(int n, int k, double p)
{
  // cumulative: P(Binom(n, p) >= k)
  if (n < 0 || k < 0)
    return 0;
  // mid p value correction
  double pb = binomial(n, k, p) / 2;
  for (int i = k + 1; i <= n; i++){
    pb += binomial(n, i, p);
  }
  return pb;
}

static bool abs_compare(int a, int b)
{
    return (std::abs(a) < std::abs(b));
}

int BETfunction::getStats()
{
	return Stats;
}

double BETfunction::getPvalue()
{
  return pvalue;
}

vector<string> BETfunction::getInteraction()
{
	return interaction_str;
}

vector<int> BETfunction::getSymmStats()
{
	return out_symmstats;
}

vector<vector<string>> BETfunction::getSymmInteraction()
{
	return out_symminter;
}

vector<string> BETfunction::getBinary()
{
  return binary_inter;
}

BETfunction:: BETfunction(NumericMatrix& X_R, int depth, bool unif, bool asymptotic)
{
	unifMargin = unif;
  asympt = asymptotic;
  X = imp(X_R);
	d = depth;
	size_t n = X.size();
	p = X[0].size();


	vector<vector<int>> c = ecdf_loc(X);
	// create a map that count for observations in the same location
	map<vector<int>, int> mapC = groupC(c);
	size_t numCount = mapC.size();

	// number of each location
	vector<int> countValues;
	for (map<vector<int>,int>::iterator it=mapC.begin(); it!=mapC.end(); ++it){
		countValues.push_back(it->second);
	}

	bid = BIDs();

	vector<string> inter = interaction_index(1);
	vector<string> inter_nb = interaction_index(0);

	// 3-dimension matrix: Xij ~ BIDs
	vector<vector<vector<int>>> matrix = CBIDs(mapC);

	// all variables go over all bids:
	vector<vector<size_t>> allidx = allComb();
	size_t total = (int)round(pow(2, p*d));


	// store symmetry statistics
	out_symmstats.resize(total, 0);

	// store interactions
	// symminter.resize(p*total);
	vector<string> temp(total);
	out_symminter.resize(p, temp);

	// store binary interactions
	binary_inter.resize(total);

	// start multiprocessing
	// loops: 2^pd
	// steps of multi-processing
	size_t remainder = total % numThread;
	size_t step = total / numThread;

	// # of assignments each thread
	vector<size_t> interval(numThread, step);
	for (size_t i = 0; i < remainder; i++){
		interval[i] += 1;
	}

	vector<size_t> thread(numThread + 1);
	for (size_t i = 1; i < numThread + 1; i++){
		thread[i] = thread[i - 1] + interval[i - 1];
	}

#ifdef _OPENMP
	omp_set_num_threads(numThread);
#pragma omp parallel for
#endif
	for (size_t th = 1; th < numThread + 1; th++){
//		double wtime = omp_get_wtime();
		for (size_t i = thread[th - 1]; i < thread[th]; i++){
		  string bi = "";
		  for (size_t v = 0; v < p; v++){
		    // store interaction
		    out_symminter[v][i] = inter[allidx[i][v]];
		    bi += inter[allidx[i][v]];
		  }
		  binary_inter[i] = bi;
			if (count(allidx[i].begin(), allidx[i].end(), 0) <= (int)(p - 2)){
				int loc0 = 1, sum = 0;
				for (size_t j = 0; j < numCount; j++){
					for (size_t k = 0; k < p; k++){
						loc0 = (~(loc0 ^ matrix[k][allidx[i][k]][j])) & 1;
					}
					sum += loc0 * countValues[j];
					loc0 = 1;
				}
				//store symmetry statistics
				out_symmstats[i] = 2 * sum - (int)n;
			}
		}
//		wtime = (omp_get_wtime() - wtime);
//		cout << "thread " + to_string(omp_get_thread_num()) + " time: " + to_string(wtime) +"\n";
	}

	// find the extreme stat
	vector<int>::iterator findMax = max_element(out_symmstats.begin(), out_symmstats.end(), abs_compare);
	size_t max_idx = distance(out_symmstats.begin(), findMax);
	// maxStat
	Stats = out_symmstats[max_idx];
	// max interaction
	interaction_str.resize(p);
	// store marginal
	int Sa = 0, Sb = 0;

	for (size_t i = 0; i < p; i++){
	  vector<string>::iterator ans = find(inter.begin(), inter.end(), out_symminter[i][max_idx]);
	  int index = distance(inter.begin(), ans);
	  interaction_str[i] = inter_nb[index];
	  if (p == 2 && i == 0){
	    // count marginal
	    for (size_t j = 0; j < numCount; j++){
	      Sa += matrix[i][index][j] * countValues[j];
	    }
	  }
	  if (p == 2 && i == 1){
	    // count marginal
	    for (size_t j = 0; j < numCount; j++){
	      Sb += matrix[i][index][j] * countValues[j];
	    }
	  }
	}

	if (p > 2 || asympt){
	  // asymptotic p value: normal N (0, n)
	  double z = abs(Stats)/sqrt((double)n);
	  pvalue = (1-0.5 * erfc(-z * sqrt(0.5))) * ((int)round(pow(2, (int)p*d)) - (int)p*((int)round(pow(2, d)) - 1) - 1) * 2;
	}else if (unifMargin || p == 1){
	  // p value: if unif of p = 1 >> binomial
	  // bonferroni
	  pvalue = ((int)round(pow(2, (int)p*d)) - (int)p*((int)round(pow(2, d)) - 1) - 1) * pbinom(n, (Stats + n) / 2, 1/2);
	}else if (p == 2){
	  // p = 2: fisher exact
	  int n11 = ((n + Stats) / 2 + Sb - (n - Sa))/2;
	  int n00 = (n + Stats) / 2 - n11;
	  int n10 = Sa - n11;
	  int n01 = Sb - n11;
	  pvalue = ((int)round(pow(2, (int)p*d)) - (int)p*((int)round(pow(2, d)) - 1) - 1) * fisherExact(n11, n01, n10, n00);
	}

	if (pvalue > 1){
	  pvalue = 1;
	}
	if (pvalue < 0){
	  pvalue = 0;
	}
}

//[[Rcpp::export]]
List symmCpp(NumericMatrix& X, int d, bool unif)
{
  BETfunction bet(X, d, unif, 1);
  size_t p = X.ncol();
  // size_t length = (int)round(pow(2, (int)p*d));

  DataFrame symm = DataFrame::create(
    Named("SymmetryStatistics") = bet.getSymmStats()
  );

  // for (size_t i = p; i > 1; i--){
  //   vector<string>::const_iterator first = bet.getSymmInteraction().end() - length*(p-i+1);
  //   vector<string>::const_iterator last = bet.getSymmInteraction().end() - length*(p-i);
  //   vector<string> temp(first, last);
  //   symm.push_front(temp, "X" + to_string(i));
  // }
  //
  // vector<string> temp1(length);
  // for (size_t i = 0; i < length; i++){
  //   temp1[i] = bet.getSymmInteraction()[i];
  // }
  // symm.push_front(temp1, "X1");

  for (size_t i = p; i > 0; i--){
    symm.push_front(bet.getSymmInteraction()[i-1], "X" + to_string(i));
  }

  symm.push_front(bet.getBinary(), "Binary Index");

  List L = List::create(Named("SymmetryStatistics") = DataFrame(symm));

  return L;

}

//[[Rcpp::export]]
List BETCpp(NumericMatrix& X, int d, bool unif, bool asymptotic)
{
  BETfunction bet(X, d, unif, asymptotic);
  int n = X.nrow();
  double zstats = abs(bet.getStats())/sqrt(n);

  DataFrame df = DataFrame::create();

  size_t p = X.ncol();
  for (size_t i = p; i > 0; i--){
    CharacterVector vi = {bet.getInteraction()[i-1]};
    df.push_front(bet.getInteraction()[i-1], "X" + to_string(i));
  }

  List L = List::create(Named("Interaction") = DataFrame(df), Named("Extreme.Asymmetry") = bet.getStats(), Named("p.value.bonf") = bet.getPvalue(), Named("z.statistic") = zstats);

  return L;

}


