/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdsize
Empirical::size() const {
  EmpiricalBuffer* p = (EmpiricalBuffer*)buffer;
  return p->ndata;
}
