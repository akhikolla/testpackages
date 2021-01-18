/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 * Must return the largest floating point values, since
 * several rank-based methods assume it.
 */
mdreal
medusa::rnan() {
  if(sizeof(mdreal) < sizeof(double)) return FLT_MAX;
  if(sizeof(mdreal) < sizeof(long double)) return DBL_MAX;
  return LDBL_MAX;
}
