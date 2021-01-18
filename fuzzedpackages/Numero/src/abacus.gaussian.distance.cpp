/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdreal
Gaussian::distance(const mdreal f, const mdreal xmu,
		   const mdreal xsigma) const {
  mdreal rlnan = medusa::rnan();
  if(f == rlnan) return rlnan;
  if(xmu == rlnan) return rlnan;
  if(xsigma == rlnan) return rlnan;

  /* Collect sampling points. */
  vector<mdreal> x;
  mdsize nloci = qloci.size();
  for(mdsize k = 0; k < nloci; k++)
    x.push_back(values[qloci[k]]);

  /* Apply transformation. */
  this->apply(x, f);

  /* Evaluate Z-scores. */
  long double d = 0.0;
  for(mdsize k = 0; k < qloci.size(); k++) {
    mdreal zk = (x[k] - xmu)/(xsigma + 1e-9);
    mdreal zdelta = (zk - zscores[k]);
    d += fabs(zdelta);
  }
  return d;
}
