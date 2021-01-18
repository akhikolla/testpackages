/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
Minimizer::Minimizer() {
  this->buffer = new MinimizerBuffer();
}

/*
 *
 */
Minimizer::Minimizer(const Minimizer& t) {
  this->buffer = new MinimizerBuffer(t.buffer);
}

/*
 *
 */
void
Minimizer::operator=(const Minimizer& t) {
  if(this == &t) return;
  MinimizerBuffer* p = (MinimizerBuffer*)buffer; delete p;
  this->buffer = new MinimizerBuffer(t.buffer);
}

/*
 *
 */
Minimizer::~Minimizer() {
  MinimizerBuffer* p = (MinimizerBuffer*)(this->buffer);
  delete p;
}
