// MIT License
// Copyright (c) 2019 Marco Pascucci

#include "linalg.h"
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

template <typename Tx, typename Ty> struct Point {
  Tx x;
  Ty y;
  Point();
  Point(Tx x, Ty y) : x(x), y(y){};
  Point(Point<Tx, Ty> &p) : x(p.x), y(p.y){};
  void operator=(Point<Tx, Ty> &p) {
    x = p.x;
    y = p.y;
  }
};

template <typename Tx, typename Ty> struct PeltResult {
  vector<unsigned int> cp;
  vector<Tx> x;
  vector<Ty> y;
  double cost;
  PeltResult(){};
  PeltResult<Tx, Ty>(vector<unsigned int> &cp, vector<Tx> &x, vector<Ty> &y,
                     double cost)
      : cp(cp), x(x), y(y), cost(cost){};
};

template <typename T1, typename T2>
double sum_of_squares(vector<T1> &y, vector<T2> &y0, size_t start, size_t end) {
  /* @brief Calculate sum of square errors between y and y0 */
  if (y.size() != y0.size()) {
    throw std::invalid_argument("a and b must have the same dimension");
  }
  if (start < 0 || end > y.size()) {
    throw std::invalid_argument("invalid start/end range");
  }
  double sum_of_squares = 0;
  for (size_t i = start; i < end; i++) {
    sum_of_squares = sum_of_squares + pow((double)(y[i] - y0[i]), 2);
  }
  return sum_of_squares;
}
template <typename T1, typename T2>
double sum_of_squares(vector<T1> &y, vector<T2> &y0) {
  return sum_of_squares(y, y0, 0, y.size());
}

template <typename T1, typename T2>
double cost_linear(vector<T1> &x, vector<T2> &y, double coeff, double inter,
                   size_t start, size_t end, bool average = false) {
  /* calculate the sum of square errors between y and x*coeff+inter */
  assert(y.size() == x.size());
  assert(start >= 0 && end <= x.size());

  double cost = 0;
  for (unsigned int i = start; i < end; i++) {
    cost = cost + pow((y[i] - (x[i] * coeff + inter)), 2);
  }
  if (average)
    cost = cost / x.size();
  return cost;
}
template <typename T1, typename T2>
double cost_linear(vector<T1> &x, vector<T2> &y, double coeff, double inter,
                   bool average = false) {
  return cost_linear(x, y, coeff, inter, 0, x.size(), average);
}

vector<unsigned int> backtrack(vector<unsigned int> cp_raw) {
  /* @brief  backtrack the changepoint vector */
  unsigned int x = cp_raw[cp_raw.size() - 1];
  vector<unsigned int> cp;
  while (x > 0) {
    cp.push_back(x);
    x = cp_raw[x];
  }
  reverse(cp.begin(), cp.end());
  return cp;
}

template <typename Tx, typename Ty>
PeltResult<Tx, Ty> op2D(vector<Tx> &x, vector<Ty> &y, double beta) {
  /* detect changepoints with PELT algorithm */
  assert(x.size() == y.size());
  unsigned int n = y.size();
  unsigned int tau;
  vector<double> q(n);
  vector<unsigned int> P(1, 0);
  vector<unsigned int> temp;
  vector<double> coeffs(n); // for result output
  vector<double> inters(n); // for result output
  vector<unsigned int> cp_raw(n, 0);
  double coeff_temp, inter_temp; // for result output

  double coeff, inter, this_cost, q_temp;
  double cp_temp = 0;
  vector<double> q_new;
  // cp_raw[0] = 0;

  for (size_t t = 1; t < n; ++t) {
    // possible changepoint @ t

    lin_reg(x, y, &coeff, &inter, 0, t + 1);
    q_temp = cost_linear(x, y, coeff, inter, 0, t + 1);
    // cp_temp = 0;
    q_new = vector<double>(t, 0);
    q_new[0] = q_temp;
    coeff_temp = coeff; // for result output
    inter_temp = inter; // for result output

    // return vector<int>(1,1);

    for (unsigned int i = 0; i < P.size(); i++) {
      tau = P[i];
      if (t - tau == 1)
        continue;
      lin_reg(x, y, &coeff, &inter, tau + 1, t + 1);
      this_cost = cost_linear(x, y, coeff, inter, tau + 1, t + 1);
      q_new[tau] = q[tau] + this_cost; // q_new
      if (q_new[tau] + beta < q_temp) {
        q_temp = q_new[tau] + beta;
        cp_temp = tau;
        coeff_temp = coeff; // for result output
        inter_temp = inter; // for result output
      }
    }
    q[t] = q_temp;
    cp_raw[t] = cp_temp;
    coeffs[t] = coeff_temp; // for result output
    inters[t] = inter_temp; // for result output

    temp.clear();

    for (size_t i = 0; i < P.size(); i++) {
      if (q_new[P[i]] < q_temp) {
        // printf("erased: %d", i);
        temp.push_back(P[i]);
      }
    }
    temp.push_back(t);
    P = temp;
    // printv(P);
  }

  // construct result structure
  vector<unsigned int> cp = backtrack(cp_raw);
  vector<unsigned int> cp1(cp);
  cp1.insert(cp1.begin(), 0);
  cp1.push_back(n - 1);
  vector<Ty> ys((cp1.size() - 1) * 2);
  vector<Tx> xs((cp1.size() - 1) * 2);

  Tx xa, xb;
  Ty ya, yb;

  for (unsigned int i = 0; i < cp1.size() - 1; i++) {
    if (i == 0) {
      xa = x[cp1[i]];
    } else {
      xa = x[cp1[i]] + 1;
    }

    xb = x[cp1[i + 1]];

    ya = coeffs[cp1[i + 1]] * xa + inters[cp1[i + 1]];
    yb = coeffs[cp1[i + 1]] * xb + inters[cp1[i + 1]];

    xs[i * 2] = xa;
    xs[i * 2 + 1] = xb;
    ys[i * 2] = ya;
    ys[i * 2 + 1] = yb;
  }

  return PeltResult<Tx, Ty>(cp, xs, ys, q[y.size() - 2]);
}
