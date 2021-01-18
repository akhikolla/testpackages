/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
Model::Model() {
  this->buffer = new ModelBuffer();
}

/*
 *
 */
Model::Model(const Topology& topo, const mdsize nsub, const mdreal eq) {
  ModelBuffer* p = new ModelBuffer();
  p->ntrain = nsub;
  p->equality = eq;
  p->topology = topo;
  this->buffer = p;
}

/*
 *
 */
Model::Model(const Model& a) {
  this->buffer = new ModelBuffer(a.buffer);
}

/*
 *
 */
void
Model::operator=(const Model& t) {
  if(this == &t) return;
  ModelBuffer* p = (ModelBuffer*)buffer; delete p;
  this->buffer = new ModelBuffer(t.buffer);
}

/*
 *
 */
Model::~Model() {
  ModelBuffer* p = (ModelBuffer*)buffer;
  delete p;
}

/*
 *
 */
mdsize
Model::order() const {
  ModelBuffer* p = (ModelBuffer*)buffer;
  return (p->codebook).order();
}

/*
 *
 */
mdsize
Model::size() const {
  ModelBuffer* p = (ModelBuffer*)buffer;
  return (p->points).size();
}

