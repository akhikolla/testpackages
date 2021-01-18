/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
Empirical::Empirical() {
  EmpiricalBuffer* p = new EmpiricalBuffer();
  this->buffer = p;
}

/*
 *
 */
Empirical::Empirical(const Empirical& t) {
  this->buffer = new EmpiricalBuffer(t.buffer);
}

/*
 *
 */
void
Empirical::operator=(const Empirical& t) {
  if(this == &t) return;
  EmpiricalBuffer* p = (EmpiricalBuffer*)buffer; delete p;
  this->buffer = new EmpiricalBuffer(t.buffer);
}

/*
 *
 */
Empirical::~Empirical() {
  EmpiricalBuffer* p = (EmpiricalBuffer*)buffer;
  delete p;
}
