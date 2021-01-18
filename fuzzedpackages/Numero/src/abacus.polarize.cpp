/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 * Convert cartesian line to polar coordinates (distance, angle).
 */
pair<mdreal, mdreal>
abacus::polarize(const mdreal x0, const mdreal y0,
		 const mdreal x, const mdreal y) {
  mdreal rlnan = medusa::rnan();
  pair<mdreal, mdreal> result(rlnan, rlnan);
  if(x0 == rlnan) return result;
  if(y0 == rlnan) return result;
  if(x == rlnan) return result;
  if(y == rlnan) return result;
  long double dx = (x - x0);
  long double dy = (y - y0);
  result.first = sqrt(dx*dx + dy*dy);
  if(dx > 0) result.second = atan(dy/(dx + 1e-10));
  else result.second = (atan(dy/(dx - 1e-10)) + M_PI);
  return result;
}
