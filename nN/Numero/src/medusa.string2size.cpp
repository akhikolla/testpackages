/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
mdsize
medusa::string2size(const string& x) {
  mdsize n = x.size();
  mdreal sznan = medusa::snan();

  /* Parse integer value. */
  long y; 
  if(sizeof(mdsize) <= sizeof(int)) y = atoi(x.c_str());
  if(sizeof(mdsize) == sizeof(long)) y = atol(x.c_str());
  if(y > 0) return y;

  /* Check preceding non-digits. */
  mdsize nsign = 0;
  mdsize nzero = 0;
  for(mdsize i = 0; i < n; i++) {
    char c = x[i];
    if(nsign > 1) return sznan;
    if(c == '+') {nsign++; continue;}
    if(isspace(c)) continue;
    if(c != '0') return sznan;
    nzero++;
  }
  if(nzero > 0) return 0;
  return sznan;
}
