/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdreal
Empirical::p(const mdreal x, const int flag) const {
  EmpiricalBuffer* p = (EmpiricalBuffer*)buffer;
  mdreal rlnan = medusa::rnan();
  long double qsmall = 0.0;
  long double qlarge = 0.0;
  if(x == rlnan) return rlnan;

  /* Count cumulative sums. */
  unordered_map<mdreal, mdreal>::iterator pos;
  unordered_map<mdreal, mdreal>& hits = p->data;
  for(pos = hits.begin(); pos != hits.end(); pos++) {
    if(pos->first <= x) qsmall += pos->second;
    if(pos->first >= x) qlarge += pos->second;
  }

  /* Normalize to unit sum. */
  double qsum = (qsmall + qlarge);
  if(qsum < 1e-20) return rlnan;
  qsmall /= qsum;
  qlarge /= qsum;

  /* Return results. */
  if(flag < 0) return qsmall;
  if(flag > 0) return qlarge;
  if(qsmall < qlarge) return 2*qsmall;
  return 2*qlarge;
}
