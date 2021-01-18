/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdreal
Gaussian::quality() const {
  
  /* Check if enough variation. */
  mdsize nuniq = 0;
  mdsize nelem = values.size();
  for(mdsize i = 1; i < nelem; i++)
    nuniq += (values[i] > values[i-1]);
  if(nuniq < (10 + sqrt(nelem))) return -1.0;

  /* Estimate quality score. */
  mdsize nloci = qloci.size();
  mdreal d = this->distance(factor, mu, sigma);
  if(d == medusa::rnan()) return 0.0;
  return nloci/(d + nloci + 1e-9);
}
