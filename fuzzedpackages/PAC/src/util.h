#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <tgmath.h>
#include <numeric>
#include <random>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "type.h"

using std::vector;
using std::string;


template<typename T>
T vecmax(const vector<T>& v) {
	typename vector<T>::const_iterator it = max_element(v.begin(), v.end());
	return *it; 
}

template<typename T>
T vecsum(const vector<T>& v) {
	T s = 0;
  for(uint i = 0; i < v.size(); i++) {
    s += v[i];
  }
	return s;
}


double reclnArea(const matrix &r);
double vecPartialSum(const vector<double>& v, const uint end);
uint whichMax(const vector<double>& vec);
int randsample(const int l, const int u, const vector<double> & p);
int randsample(const int l, const int u, const vector<double> & p);

uint num_unique(vector<int> vec);
int nwithin(const matrix &range, const matrix *data);
bool within(const matrix &range, const vector<double> & data);

#endif /* UTIL_H */
