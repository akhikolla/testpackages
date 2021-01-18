/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 * Convert integers into human readable form.
 */
string
medusa::long2text(const long value) {
  char buffer[32];

  /* Small value. */
  long magn = labs(value);
  if(magn < 1000) {
    sprintf(buffer, "%d", (int)value);
    return string(buffer);
  }

  /* Collect thousand-sets. */
  vector<int> segments;
  unsigned long k = magn;
  while(k > 0) {
    segments.push_back(k%1000);
    k /= 1000;
  }

  /* Finish text. */
  string s; int sign = 1;
  if(value < 0) sign = -1;
  for(int i = (int)(segments.size() - 1); i >= 0; i--) {
    if(s.size() < 1) sprintf(buffer, "%d", sign*segments[i]);
    else sprintf(buffer, ",%03d", segments[i]);
    s += string(buffer);
  }
  return s;
}
