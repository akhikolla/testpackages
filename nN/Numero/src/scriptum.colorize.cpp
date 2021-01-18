/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 * Return red, green and blue components.
 */
vector<Color>
scriptum::colorize(const vector<mdreal>& y, const mdreal scale,
		   const std::string& name) {
  mdsize nelem = y.size();
  mdreal rlnan = medusa::rnan();
  
  /* Calculate scaled statistics. */
  mdreal center = statistic(y, "center");
  mdreal sigma = statistic(y, "sd");

  /* Null color. */
  vector<Color> colors(nelem);
  if(sigma == rlnan) return colors;

  /* Determine pseudo-coloring. */
  if(sigma < 1e-10) sigma = 1e-10;
  for(mdsize i = 0; i < nelem; i++) {
    if(y[i] == rlnan) continue;
    double value = scale*(y[i] - center)/sigma;
    colors[i] = colormap((0.25*value + 0.5), name);
  }
  return colors;
}
