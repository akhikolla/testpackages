
#include "Bounds.h"

#include <Rversion.h>
#if defined(R_VERSION) && R_VERSION >= R_Version(3, 3, 0)
  /* nothing */
#else
  #define NO_C_HEADERS true // disables the including of stdlib.h, stdio.h, limits.h, math.h by R.h
  #include <cstdlib>        // manually loading of cpp versions of disabled headers
  #include <cstdio>
  #include <climits>
  #include <cmath>
  #include <cstddef>
  using std::size_t;
#endif

#include <R.h>
#include <Rmath.h>

/***************
* class LUBound
* specifies lower and upper bound
* Thomas Hotz, 2007-2011
***************/

/*************
* constructor given lower and upper bound
****************/
LUBound::LUBound(double lb, double ub) : lower(lb), upper(ub) {}

/*************
* constructor with default bounds set to -+Inf, i.e. always fulfilled
****************/
LUBound::LUBound() : lower(R_NegInf), upper(R_PosInf) {}

/*************
* add
* add more constraints
****************/
void LUBound::add(double lb, double ub) {
  lower = fmax2(lower, lb);
  upper = fmin2(upper, ub);
}
void LUBound::add(LUBound b) {
  add(b.lower, b.upper);
}

/*************
* feasible
* check feasibility of constraints
****************/
bool LUBound::feasible() {
  return lower <= upper;
}

/***************
 * class Bounds
 * stores lower and upper bounds over intervals
 * Thomas Hotz, 2007-2011
 ***************/

/*************
* constructor
****************/
Bounds::Bounds(unsigned int n, int* xli, unsigned int ni, int* xri, double* xlb, double * xub) : N(n), li(xli), Ni(ni), ri(xri), lb(xlb), ub(xub) {
  // ensure there are bounds
  if(Ni < 1) error("no bounds specified!");

  // allocate arrays
  nexti = (int*) R_alloc(N, sizeof(int));
  cri = (int*) R_alloc(N, sizeof(int));
  cb = (LUBound*) R_alloc(N, sizeof(LUBound));

  // initialize ci, cri and cb, and check whether bounds can be fulfilled; corresponds to k = 0
  for(unsigned int i = 0; i < N; i++) {
    // initialize bound to be infinite
    cb[i] = LUBound();
    // current right index equals left index
    cri[i] = i;
    // go through all intervals [i,i]
//     Rprintf("i = %d, nexti[i] = %d, cri[i] = %d\n", i, nexti[i], cri[i]);
    for(nexti[i] = li[i]; nexti[i] != R_NaInt && nexti[i] < (int) Ni && ri[nexti[i]] == (int) i; nexti[i]++) {
//       Rprintf("  i = %d, nexti[i] = %d, ri[nexti[i]] = %d\n", i, nexti[i], ri[nexti[i]]);
      // check whether we've gone too far
      if(i < N - 1 && li[i+1] != R_NaInt && nexti[i] >= li[i+1]) {
	nexti[i] = R_NaInt; // indicate that there are no more constraints for this left index
// 	Rprintf("  break: i = %d, nexti[i] = %d, li[i+1] = %d\n", i, nexti[i], li[i+1]);
	break;
      }
      cb[i].add(lb[nexti[i]], ub[nexti[i]]); // add constraint
//       Rprintf("  i = %d, cb[i].lower = %4.3e, cb[i].upper = %4.3e\n", i, cb[i].lower, cb[i].upper);
    }
    // check whether we've gone too far
    if(nexti[i] >= (int) Ni) {
      nexti[i] = R_NaInt; // indicate that there are no more constraints for this left index
    }
//     // check well-behaviour of indices
//     if(nexti[i] != R_NaInt && nexti[i] >= Ni) error("index %d of interval with left index %d is >= number of intervals %d!", nexti[i], i, Ni);
    // check feasibility
    #ifdef DEBUGbounded
    Rprintf("i = %d, cb[i].lower = %4.3e, cb[i].upper = %4.3e\n", i, cb[i].lower, cb[i].upper);
    #endif
    if(!cb[i].feasible()) error("bounds not feasible at index %d!", i);
  }
}

/*************
* current
* return current bound for interval [l,r], assuming it has been computed for [l,r-1] and [l+1,r] before
****************/
LUBound Bounds::current(unsigned int l, unsigned int r) {
  // check whether these indices may be asked for
  if(l >= N || r >= N  || r < l) error ("indices must fulfill l %d <= r %d < N %d", l, r, N);
  if((int) r < cri[l]) error("for l %d we are already at cri %d, i.e. beyond r %d", l, cri[l], r);
  if((int) r > cri[l] + 1) error("for l %d we are at cri %d, i.e. r %d is too far", l, cri[l], r);
  // check whether we have already computed this one
  if(cri[l] == (int) r) return cb[l];
  // add already computed bound on [li+1,ri]
  if(l < N - 1 && cri[l + 1] != (int) r) {
    error("bound for l + 1 = %d and r = %d needs to be available, but is at cri %d!", l + 1, r, cri[l + 1]);
  } else {
    cb[l].add(cb[l + 1]); // add constraint
  }
  // add all intervals [l,r] in the list
  for(; nexti[l] != R_NaInt && nexti[l] < (int) Ni && ri[nexti[l]] == (int) r; nexti[l]++) {
//       Rprintf("  l = %d, nexti[l] = %d, ri[nexti[l]] = %d\n", l, nexti[l], ri[nexti[l]]);
      // check whether we've gone too far
      if(l < N - 1 && li[l+1] != R_NaInt && nexti[l] >= li[l+1]) {
	nexti[l] = R_NaInt; // indicate that there are no more constraints for this left index
// 	Rprintf("  break: l = %d, nexti[l] = %d, li[l+1] = %d\n", l, nexti[l], li[l+1]);
	break;
      }
      cb[l].add(lb[nexti[l]], ub[nexti[l]]); // add constraint
//       Rprintf("  l = %d, cb[l].lower = %4.3e, cb[l].upper = %4.3e\n", l, cb[l].lower, cb[l].upper);
  }
  // update current right index
  cri[l] = r;
  // return bound on this interval
  #ifdef DEBUGbounded
  Rprintf("l = %d, r = %d, cb[l].lower = %4.3e, cb[l].upper = %4.3e\n", l, r, cb[l].lower, cb[l].upper);
  #endif
  return cb[l];
}
