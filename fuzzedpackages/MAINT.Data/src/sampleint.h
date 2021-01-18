#ifndef SAMPLEINTH
#define SAMPLEINTH

//#include <Rcpp.h>
#include "RcppArmadillo.h"
#include <vector>

using namespace Rcpp;
using namespace std;

template<class T>
struct Comp{
  Comp(const vector<T>& v ) : _v(v) {}
  bool operator ()(int a, int b) { return _v[a] > _v[b]; }
  const vector<T>& _v;
};

template<class T>
void highestval(const unsigned n,const int k,const vector<T>& val,
  vector<int>& res,const bool C0cnv)
{
  static vector<int> vx;
  if (vx.size()!=n) vx.resize(n);
  for(unsigned i=0; i<n; ++i) vx[i] = i;
  partial_sort(vx.begin(), vx.begin()+k, vx.end(), Comp<T>(val));
  for(int i=0; i<k; ++i) {
    if (C0cnv) res[i] = vx[i];
    else res[i] = vx[i] + 1;
  }
  return;
};

void sampleint(const int n,const int size,vector<int>& res,const bool C0cnv);

#endif
