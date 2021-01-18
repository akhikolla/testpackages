/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
mdreal
medusa::rlim() {
  if(sizeof(mdreal) < sizeof(double)) return FLT_MIN;
  if(sizeof(mdreal) < sizeof(long double)) return DBL_MIN;
  return LDBL_MIN;
}
