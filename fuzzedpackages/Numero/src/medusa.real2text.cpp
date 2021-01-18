/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
string
medusa::real2text(const mdreal x) {
  mdreal ipart = 0.0;
  char buf[32];
  if(x == medusa::rnan()) return "nan";
  if(x == 0.0) return "0";
  mdreal amp = fabs(x);
  if(amp > 1e24) {return real2string(x);}
  if(amp > 4.99e9) {sprintf(buf, "%+.0fe9", x/1e9); return string(buf);};
  if(amp > 4.99e6) {sprintf(buf, "%+.0fe6", x/1e6); return string(buf);};
  if(amp > 4999.5) {sprintf(buf, "%+.0fe3", x/1000.0); return string(buf);};
  if((amp > 14.9) || (modf(amp, &ipart) == 0.0)) {
    sprintf(buf, "%+.0f", x);
    return string(buf);
  };
  if(amp > 4.99) {sprintf(buf, "%+.1f", x); return string(buf);};
  if(amp > 0.499) {sprintf(buf, "%+.2f", x); return string(buf);};
  if(amp > 0.0499) {sprintf(buf, "%+.3f", x); return string(buf);};
  if(amp > 0.00499) {sprintf(buf, "%+.4f", x); return string(buf);};
  sprintf(buf, "%+.2e", x);
  return string(buf);
}
