/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
string
medusa::long2string(const long value) {
  char buf[128];
  sprintf(buf, "%ld", value);
  return string(buf);
}
