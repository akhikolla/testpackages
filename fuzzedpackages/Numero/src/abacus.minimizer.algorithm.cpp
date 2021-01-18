/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
void
Minimizer::algorithm(const mdsize nbins, const mdreal toler) {
  MinimizerBuffer* p = (MinimizerBuffer*)buffer;
  mdreal rlnan = medusa::rnan();
  p->epsilon = toler;
  if(p->epsilon < 0.0) p->epsilon = 0.0;
  if(p->epsilon == rlnan) p->epsilon = 0.0;
  p->npoints = nbins;
}

/*
 *
 */
std::pair<mdsize, mdreal>
Minimizer::algorithm() const {
  MinimizerBuffer* p = (MinimizerBuffer*)buffer;
  return std::pair<mdsize, mdreal>(p->npoints, p->epsilon);
}
