#ifndef SPECIAL_FUNCTIONS_H
#define SPECIAL_FUNCTIONS_H

#include <algorithm>

#include <Rmath.h>

namespace {
  inline double safe_product(const double x,
			     const double y) {
    const double xy = x * y;
    if ( (x==0) && ( (y == R_PosInf) || (y == R_NegInf)) ) {
      return 0;
    } else {
      return std::max(xy, -1.0);
    }
  }
  
  inline
  double dexprl(const double x) {
    if (R_FINITE(x)) {
      /* actual value */
      if (x != 0) {
	return expm1(x) / x;
      } else {
	return 1;
      }
    }
    else if (ISNA(x)) {
      return NA_REAL;
    }
    else if (ISNAN(x)) {
      return R_NaN;
    }
    else if (x == R_PosInf) {
      return R_PosInf;
    }
    else if (x == R_NegInf) {
      return 0;
    } else {
      /* can't happen */
      return NA_REAL;
    }
  }
  
  inline
  double log1prel(const double x) {
    if (R_FINITE(x)) {
      /* actual value */
      if (x != 0) {
	return log1p(x) / x;
      } else {
	return 1;
      }
    }
    else if (ISNA(x)) {
      return NA_REAL;
    }
    else if (ISNAN(x)) {
      return R_NaN;
    }
    else if (x == R_PosInf) {
      return 0;
    }
    else if (x == R_NegInf) {
      return R_NaN;
    } else {
      /* can't happen */
      return NA_REAL;
    }
  }
  
  /*
    if someone could tell me why this isn't a part of base R
    I'd be insanely grateful
  */
  
  inline
  double texmex_log1mexp(const double x) {
    /* accurately compute log(1-exp(x)) */
    if (R_FINITE(x)) {
      if (x > -M_LN2) {
	return log(-expm1(x));
      } else {
	return log1p(-exp(x));
      }
    }
    else if (ISNA(x)) {
      return NA_REAL;
    }
    else if (ISNAN(x)) {
      return R_NaN;
    }
    else if (x == R_PosInf) {
      return R_NaN;
    }
    else if (x == R_NegInf) {
      return 0;
    } else {
      /* can't happen */
      return NA_REAL;
    }
  }
}

#endif
