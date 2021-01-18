// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
using std::vector;

#include <numeric> // accumulate
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/tools/roots.hpp>
#include <iostream> // cout
//#include <algorithm>
//#include <random>
//#include <functional>
//#include <boost/array.hpp>
//#include <boost/math/distributions/normal.hpp>
//#include <math.h>       /* round, floor, ceil, trunc */
//using namespace boost::math::tools;           // For bracket_and_solve_root.
using namespace boost::math::tools;           // For bracket_and_solve_root.
std::vector<double> cmix_g;
int m_g;
double alp;




bool isZero(double i)
{
  return i == 0;
}

// function to remove 0 and to sort the vector in descending order //
void remo(std::vector<double> &v1) {
  std::vector<double>::iterator newIter = remove_if(v1.begin(), v1.end(), isZero);
  v1.resize(newIter - v1.begin());
}

// function to weed extremely small eigenvalues

std::vector<double> weed(std::vector<double> v2, double acc = 50.0) { // v here is without zeros
  sort(v2.rbegin(), v2.rend()); // sort elements in decreasing order
  int mw = v2.size();
  int qw;
  int rw;
  double sw;
  double mean;
  while (mw>2 && v2[0] > 0 && v2[mw - 1] / v2[0] < (1.0 / acc)) {
    qw = mw - 1;
    rw = mw - 2;
    v2[qw - 1] += v2[mw - 1];
    while (rw > 1 && v2[rw - 1] < v2[qw - 1]) {
      sw = 0;
      for (int i = rw - 1; i <= qw - 1; ++i) {
        sw += v2[i];
      }
      mean = sw / (qw - rw + 1);
      for (int i = rw - 1; i <= qw - 1; ++i) {
        v2[i] = mean;
      }
      rw -= 1;
      
    }
    mw = qw;
  }
  v2.resize(mw);
  if (v2.size() == 2 && v2[1] / v2[0] < 1e-3) {
    v2.resize(1);
  }
  return v2;
}


// get the eigenvalue vector for a specific level in between uppperbound and lowerbound
// [[Rcpp::export]]
std::vector<double> getL(std::vector<double> ub, std::vector<double> lb, double level) {
  // lu and lb here are original eigenvalues with the same length (without removing zeros but they are both in descending order)
  int sl = lb.size();
  // the difference between upper bound and lower bound
  std::vector<double> lam1(sl);
  transform(ub.begin(), ub.end(), lb.begin(), lam1.begin(), std::minus<double>());
  // cumulative sum of the diff, with 1st value 0
  std::vector<double> cums(sl);
  cums[0] = 0;
  partial_sum(lam1.begin(), lam1.end() - 1, cums.begin() + 1);
  // the maximal diff 
  std::vector<double> lamdiff2(sl);
  double levels = accumulate(lb.begin(), lb.end(), 0.0);
  double sumd = level - levels;
  double pdiff;
  for (int i = 0; i < sl; ++i) {
    pdiff = sumd - cums[i];
    if (pdiff < 0) {
      pdiff = 0;
    }
    lamdiff2[i] = pdiff;
  }
  
  std::vector<double> mlam;
  transform(lam1.begin(), lam1.end(), lamdiff2.begin(),
            std::back_inserter(mlam), [](double a, double b) {return std::min(a, b); });
  
  std::vector<double> flam(sl);
  transform(lb.begin(), lb.end(), mlam.begin(), flam.begin(), std::plus<double>());
  return flam;
}

// function to calculate the initial beta
double ruben(std::vector<double> v3) { // v <- zero removed and weeded
  double mx = *std::max_element(v3.begin(), v3.end());
  double mn = *std::min_element(v3.begin(), v3.end());
  return 2 * mn * mx / (mn + mx);
}

// function to get c mixtures
std::vector<double> getC(std::vector<double> v, double beta, double eps = 1e-10) {// v <- zero removed and weeded
  remo(v);
  v=weed(v);
  sort(v.begin(), v.end());
  int m = v.size();
  int ix = 0;
  std::vector<double> d;
  std::vector<double> c;
  // calculate the first element of vector c
  double prod = 1;
  for (int m1 = 0; m1 < m; ++m1) {
    prod = prod * sqrt(beta / v[m1]);
  }
  c.push_back(prod);
  if (c[0] == 0 || c[0]>1)
    return{};
  double restc = 1 - c[0];
  
  std::vector<double> dbase;
  for (int m2 = 0; m2 < m; ++m2) {
    dbase.push_back(1 - beta / v[m2]);
  }
  
  bool ready = false;
  double sum3;
  double sum4;
  while (!ready) {
    sum3 = 0;
    for (int m3 = 0; m3 < m; ++m3) {
      sum3 += 0.5*pow(dbase[m3], ix + 1);
    }
    d.push_back(sum3);
    
    sum4 = 0;
    for (int i = 0; i < (ix + 1); ++i) {
      //Lam5.push_back(c[i] * d[ix - i] / (ix + 1));
      sum4 += c[i] * d[ix - i] / (ix + 1);
    }
    c.push_back(sum4);
    
    if (restc > std::numeric_limits<double>::epsilon() * 50)
      restc = restc - c[ix + 1];
    else
      restc = restc * v[0] / v[m - 1];
    
    ready = (restc < eps) || (c[ix + 1]<0);
    ix = ix + 1;
  }
  
  if (restc>0)
    c.push_back(restc);
  if (std::all_of(c.begin(), c.end(), [](double i) { return i >= 0; }))
    return c;
  else
    return{};
  
}

//functions to calculate p-value
std::vector<double> getCB(std::vector<double> lam, int length = 50) {
  remo(lam);
  lam = weed(lam);
  double minlam = *std::min_element(lam.begin(), lam.end());
  double betstart = ruben(lam);
  double dec = (betstart - minlam) / (length - 1);
  std::vector<double> cbmix = getC(lam, betstart);
  while (cbmix.size() == 0) {
    betstart -= dec;
    cbmix = getC(lam, betstart);
  }
  cbmix.push_back(betstart);
  return cbmix;
}


//functions to calculate p-value
//[[Rcpp::export]]
double pv(double x, std::vector<double> lam) {//v = original eigenvalues
  // first remove zeros and weed lams
  remo(lam);
  lam = weed(lam);
  int mp = lam.size();
  double p = 0.0;
  if (mp == 1) {
    if (lam[0] != 0) {
      boost::math::chi_squared dist1(1);
      p = 1 - cdf(dist1, x / lam[0]);
    }
    else {
      p = 1.0;
    }
  }
  else {
    std::vector<double> cmix = getCB(lam);
    int mc = cmix.size();
    for (int ci = 1; ci < mc; ++ci) {
      boost::math::chi_squared dist2(mp + 2 * (ci - 1));
      p += cmix[ci - 1] * (1 - cdf(dist2, x / cmix[mc-1]));
    }
  }
  return p;
}

////function to calculate the critical value by bracketting method
//[[Rcpp::export]]
double criticalvalue(std::vector<double> lam, double alpha=0.05) {// original v
  remo(lam);
  lam = weed(lam);
  m_g = lam.size(); 
  alp = alpha;
  double ctv = 0.0;
  if (m_g == 1) {
    if (lam[0] != 0) {
      boost::math::chi_squared dist3(1);
      ctv = quantile(dist3, 0.95) * lam[0];
    }
    else {
      ctv = 0.0;
    }
  }
  else {
    
    std::vector<double> cmixlocal = getCB(lam); 
    cmix_g = cmixlocal;
    
    struct TerminationCondition {
      bool operator() (double min_global, double max_global) {
        return abs(min_global - max_global) <= 1.0e-10;
      }
    };
    
    struct Findroot {
      double operator() (double x) {
        double sumcchi = 0;
        int cm = cmix_g.size();
        for (int kc = 1; kc < cm; ++kc) {
          boost::math::chi_squared dist4(m_g + 2 * (kc - 1));
          sumcchi += cmix_g[kc - 1] * (1 - cdf(dist4, x / cmix_g[cm - 1]));
        }
        return sumcchi - alp;// Replace with your function
      }
    };
    
    double guess = accumulate(lam.begin(), lam.end(), 0.0);
    double factor = 2.0;
    bool is_rising = false;
    const boost::uintmax_t maxit = 20;            // Limit to maximum iterations.
    boost::uintmax_t it = maxit;
    
    std::pair<double, double> r = bracket_and_solve_root(Findroot(), guess, factor, is_rising, TerminationCondition(), it);
    ctv = r.first + (r.second - r.first) / 2;
  }
  cmix_g.clear();
  return ctv;
  
}
