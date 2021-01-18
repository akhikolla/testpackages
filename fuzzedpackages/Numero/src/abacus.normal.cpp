/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
Normal::Normal() {
  Approximation* p = new Approximation();
  this->buffer = p;
}

/*
 *
 */
Normal::Normal(const Normal& t) {
  this->buffer = new Approximation(t.buffer);
}

/*
 *
 */
void
Normal::operator=(const Normal& t) {
  if(this == &t) return;
  Approximation* p = (Approximation*)buffer; delete p;
  this->buffer = new Approximation(t.buffer);
}

/*
 *
 */
Normal::~Normal() {
  Approximation* p = (Approximation*)buffer;
  delete p;
}

/*
 *
 */
bool
Normal::configure(const vector<mdreal>& prm) {
  Approximation* p = (Approximation*)buffer;
  return p->parameters(prm);
}

/*
 *
 */
vector<mdreal>
Normal::parameters() const {
  Approximation* p = (Approximation*)buffer;
  return p->parameters();
}

/*
 *
 */
mdreal
Normal::z(const mdreal x) const {
  Approximation* p = (Approximation*)buffer;
  return p->transform(x);
}

/*
 *
 */
void
Normal::z(vector<mdreal>& xdata) const {
  Approximation* p = (Approximation*)buffer;
  for(mdsize i = 0; i < xdata.size(); i++)
    xdata[i] = p->transform(xdata[i]);
}
