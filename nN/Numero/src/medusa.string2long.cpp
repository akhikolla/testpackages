/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
long
medusa::string2long(const string& x) {
  mdsize n = x.size();
  mdreal lnnan = medusa::lnan();

  /* Parse integer value. */
  long y = atol(x.c_str());
  if(y > 0) return y;

  /* Check preceding non-digits. */
  mdsize nsign = 0;
  mdsize nzero = 0;
  for(mdsize i = 0; i < n; i++) {
    char c = x[i];
    if(nsign > 1) return lnnan;
    if(c == '+') {nsign++; continue;}
    if(isspace(c)) continue;
    if(c != '0') return lnnan;
    nzero++;
  }
  if(nzero > 0) return 0;
  return lnnan;
}
