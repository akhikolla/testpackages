/* fisher.cpp - a module of OPUS Miner providing fisherTest, a function to calculate the Fisher exact test and log_combin a function to calculate the log of the number of combinations of k items selected from n.
** Copyright (C) 2012 Geoffrey I Webb
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>
#include <cmath>
#include "fisher.h"

// return the log of the factorial of n
static double logfact(const int n) {
  static std::vector<double> lf;

  int i;

  for (i = lf.size(); i <= n; i++) {
    if (i == 0) lf.push_back(0);
    else lf.push_back(lf[i-1] + std::log(static_cast<double>(i)));
  }

  return lf[n];
}

double log_combin(const int n, const int k) {
  return logfact(n) - logfact(k) - logfact(n-k);
}

// return the p value for a one tailed fisher exact test for the probability of obtaining d or more in a contingency table where the marginal frequencies are invariant
double fisherTest(int a, int b, int c, int d) {
  double p = 0; // cumulative value of p

  // will loop until b or c is 0 - as the values are interchangeable, make c the lesser value and test only for when it reaches 0
  if (b < c) {
    const int t = b;
    b = c;
    c = t;
  }

  const double invariant = -logfact(a+b+c+d)+logfact(a+b)+logfact(c+d)+logfact(a+c)+logfact(b+d);

  do {
    p += std::exp(invariant-logfact(a)-logfact(b)-logfact(c)-logfact(d));
    a++;
    b--;
    c--;
    d++;
  } while (c >= 0);

  return p;
}

