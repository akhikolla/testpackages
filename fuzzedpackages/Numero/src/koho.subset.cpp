/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

typedef map<mdreal, vector<Point*> > ContentMap;

/*
 *
 */
Subset::Subset() {
  this->label = medusa::snan();
  this->capacity = 0;
  this->occupancy = 0;
}

/*
 *
 */
Subset::~Subset() {}

/*
 *
 */
void
Subset::clear() {
  this->occupancy = 0;
  (this->contents).clear();
}

/*
 *
 */
void
Subset::configure(const mdsize rank, const mdsize cap) {
  this->label = rank;
  this->capacity = cap;
  this->clear();
}

/*
 *
 */
mdsize
Subset::size() const {
  return occupancy;
}
