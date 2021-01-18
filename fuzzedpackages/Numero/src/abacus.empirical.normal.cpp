/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
Normal
Empirical::normal() const {
  EmpiricalBuffer* p = (EmpiricalBuffer*)buffer;

  /* Make sure parameters have been estimated. */
  vector<mdreal> values;
  vector<mdreal> weights;
  p->contents(values, weights);
  (p->approx).fit(values, weights);

  /* Return parameter vector. */
  vector<mdreal> prm = (p->approx).parameters();
  Normal output; output.configure(prm);
  return output;
}
