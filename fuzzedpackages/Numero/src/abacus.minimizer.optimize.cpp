/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdreal
Minimizer::optimize(Minimizer& func) {
  mdreal rlnan = medusa::rnan();

  /* Check parameters. */
  pair<mdsize, mdreal> param = func.algorithm();
  pair<mdreal, mdreal> lims = func.space();
  if(param.first < 2) return rlnan;
  if(param.second == rlnan) return rlnan;
  if(lims.first == rlnan) return rlnan;
  if(lims.second == rlnan) return rlnan;

  /* Set starting point. */
  mdreal xmin = lims.first;
  mdreal xmax = lims.second;
  mdreal xrange = (xmax - xmin);
  mdreal xnadir = rlnan;
  mdreal ynadir = rlnan;

  /* Divide and conquer. */
  for(mdsize t = 0; t < 100; t++) {

    /* Find the point that produces the minimum value. */
    mdreal delta = (xmax - xmin)/(param.first);
    for(mdsize i = 0; i <= param.first; i++) {
      mdreal x = (xmin + i*delta);
      if(ynadir == rlnan) {
	ynadir = func.value(x);
	xnadir = x;
	continue;
      }
      if(x == xnadir) continue;
      mdreal y = func.value(x);
      if(y == rlnan) continue;
      if(y >= ynadir) continue;
      ynadir = y;
      xnadir = x;
    }
    
    /* Update search space. */
    xmin = (xnadir - delta);
    xmax = (xnadir + delta);
    if(xmin < lims.first) xmin = lims.first;
    if(xmax > lims.second) xmax = lims.second;
    
    /* Check tolerance. */
    mdreal ratio = 2*(xmax - xmin)/xrange;
    if(ratio <= param.second) return xnadir;
  }
  worry("Optimization failed.", __FILE__);
  return rlnan;
}
