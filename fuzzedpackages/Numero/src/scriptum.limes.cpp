/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
Limes::Limes() {
  this->alpha = medusa::rnan();
  this->omega = medusa::rnan();
}

/*
 *
 */
Limes::~Limes() {}

/*
 *
 */
bool
Limes::update(const mdreal x) {
  mdreal rlnan = medusa::rnan();
  if(x == rlnan) return false;
  if(x < MINCOORD_scriptum) return false;
  if(x > MAXCOORD_scriptum) return false;
  if(alpha == rlnan) alpha = x;
  if(omega == rlnan) omega = x;
  if(x < alpha) this->alpha = x;
  if(x > omega) this->omega = x;
  return true;
}

/*
 *
 */
bool
Limes::update(const mdreal x, const Style& sty) {
  mdreal rlnan = medusa::rnan();
  if(x == rlnan) return false;

  /* Check padding. */
  mdreal pad = sty.padding;
  if(pad < 0.0) pad = 0.0;

  /* Apply padding. */
  mdreal xmin = (x - pad);
  mdreal xmax = (x + pad);

  /* Update limits. */
  if(xmin < MINCOORD_scriptum) return false;
  if(xmax > MAXCOORD_scriptum) return false;
  if(alpha == rlnan) this->alpha = xmin;
  if(omega == rlnan) this->omega = xmax;
  if(xmin < alpha) this->alpha = xmin;
  if(xmax > omega) this->omega = xmax;
  return true;
}

/*
 *
 */
bool
Limes::update(const vector<mdreal>& x, const Style& sty) {
  mdreal rlnan = medusa::rnan();

  /* Find extreme values. */
  mdreal xmin = statistic(x, "min");
  mdreal xmax = statistic(x, "max");
  if(xmin == rlnan) return false;
  if(xmax == rlnan) return false;

  /* Apply padding. */
  mdreal pad = sty.padding;
  if(pad < 0.0) pad = 0.0;
  xmin -= pad;
  xmax += pad;

  /* Update limits. */
  if(xmin < MINCOORD_scriptum) return false;
  if(xmax > MAXCOORD_scriptum) return false;
  if(alpha == rlnan) this->alpha = xmin;
  if(omega == rlnan) this->omega = xmax;
  if(xmin < alpha) this->alpha = xmin;
  if(xmax > omega) this->omega = xmax;
  return true;
}
