#include <Rcpp.h>
using namespace Rcpp;

typedef ListOf<NumericVector> dblList; 
typedef std::vector<NumericVector::iterator> dblIteratorList;

// show() for dblIteratorList & std::vector<double>
void show(dblIteratorList x){
  for (size_t i=0; i<x.size(); ++i){
    NumericVector::iterator xi = x[i];
    Rcout << *xi << " ";
  }
  Rcout << std::endl;
}
void show(std::vector<double> x){
  for (size_t i=0; i<x.size(); ++i){
    Rcout << x[i] << " ";
  }
  Rcout << std::endl;
}

// getStarts and getEnds for dblIteratorList
dblIteratorList getStarts(dblList x){
  dblIteratorList out;
  int l = x.size();
  for (int i=0; i<l; ++i){
    out.push_back(x[i].begin());
  }
  return out;
}
dblIteratorList getEnds(dblList x){
  dblIteratorList out;
  int l = x.size();
  for (int i=0; i<l; ++i){
    out.push_back(x[i].end());
  }
  return out;
}

// increase for dblIteratorList
void increase(dblIteratorList& ii, int& changed,
              const dblIteratorList& _starts_, 
              const dblIteratorList& _ends_, 
              int pos = -1){
  if (pos == -1) pos = ii.size() - 1;
  ++ii[pos];
  changed = pos;
  if (ii[pos] == _ends_[pos]){
    if (pos > 0){
      ii[pos] = _starts_[pos];
    increase(ii, changed, _starts_, _ends_, pos = pos-1);
    }
  }
}

bool finished(dblIteratorList ii, 
              const dblIteratorList _ends_){
  return (ii.front() == _ends_.front());
};

// Test function
//// [[Rcpp::export]]
// void test(dblList x){
//   dblIteratorList _starts_ = getStarts(x);
//   dblIteratorList _ends_ = getEnds(x);
//   dblIteratorList ii = _starts_;
//   int changed = 0;
//   do {
//     for (int ch=changed; ch < x.size(); ++ch){
//       Rcout << "      > Inserting "<< *ii[ch] <<
//         " at position " << ch << std::endl;
//     }
//     show(ii);
//     // more actions here if required...
//     increase(ii, changed, _starts_, _ends_);
//   } while(!finished(ii, _ends_));
// }


/* R
test(list(1:2, 5:9))
test(list(4:2 * 2, 4:8, c(55, 44, 66)))
test(list(1, 2, 3:4, 5, 6:7))
*/

// double sum(dblIteratorList x){
//   double out = 0;
//   for (unsigned int i=0; i<x.size(); ++i){
//     out += *x[i];
//   }
//   return out;
// }

// Aux fun dblList -> easier using std::accumulate!
int ll(dblList x){
  int out = 1;
  for (int i=0; i<x.size(); ++i){
    out *= x[i].size();
  }
  return out;
}

// iteration using dblIteratorList
// [[Rcpp::export]]
NumericMatrix C_iterate(dblList dx, dblList dminxy, double Sx_base, double Sy,
                        bool verbose = false){
  // initialize
  dblIteratorList dx_starts = getStarts(dx), dx_ends = getEnds(dx),
    dminxy_starts = getStarts(dminxy), dminxy_ends = getEnds(dminxy);
  dblIteratorList dx_it = dx_starts, dminxy_it = dminxy_starts;
  // Define output matrix
  NumericMatrix out(ll(dx), 2);
  int count = 0, changed = 0;
  size_t n = dx.size();
  std::vector<double> cum_x(n), cum_minxy(n);
  do {
    for (size_t ch = changed; ch < n; ++ch){
      if (ch == 0){
        cum_x.front() = *dx_it.front(); cum_minxy.front() = *dminxy_it.front();
      } else {
        cum_x[ch] = cum_x[ch-1] + *dx_it[ch]; 
        cum_minxy[ch] = cum_minxy[ch-1] + *dminxy_it[ch];
      }
    }
    double numerator = Sx_base + cum_minxy.back();
    double con = numerator / (Sx_base + cum_x.back());
    double cov = numerator / Sy;
    out(count, 0) = con; out(count, 1) = cov;
    count++;
    if (verbose){ 
      Rcout << "cum_x:     "; show(cum_x);
      Rcout << "cum_minxy: "; show(cum_minxy);
      Rcout << "dx:      "; show(dx_it);
      Rcout << "dminxy:  "; show(dminxy_it);
      Rcout << "numerator=" << numerator << std::endl;
      Rcout << "con=" << con << ", cov=" << cov << std::endl;
    }
    increase(dx_it, changed, dx_starts, dx_ends);
    increase(dminxy_it, changed, dminxy_starts, dminxy_ends);
  } while(!finished(dx_it, dx_ends));
  return out;
}


