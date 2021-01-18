/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

#define MIN_SIZE_TRIGGER 8
#define SPARSE_RATIO 0.1
#define FULL_RATIO 0.2

/*
 *
 */
Array::Array() {
  this->rlnan = medusa::rnan();
  this->ndata = 0;
  this->nelem = 0;
}

/*
 *
 */
Array::~Array() {}

/*
 *
 */
mdreal
Array::operator[](const mdsize rank) const {
  if(rank >= nelem) return rlnan;
  if(full.size() > 0) return full[rank];
  map<mdsize, mdreal>::const_iterator it;
  if((it = sparse.find(rank)) == sparse.end()) return rlnan;
  return it->second;
}

/*
 *
 */
void
Array::elements(vector<Element>& elem, const mdsize rowrank) const {
  Element e;
  for(mdsize j = 0; j < full.size(); j++) {
    if((e.value = full[j]) == rlnan) continue;
    e.row = rowrank;
    e.column = j;
    elem.push_back(e);
  }
  for(map<mdsize, mdreal>::const_iterator it = sparse.begin();
      it != sparse.end(); it++) {
    e.row = rowrank;
    e.column = it->first;
    e.value = it->second;
    elem.push_back(e);
  }
}

/*
 *
 */
mdsize
Array::length() {
  return this->optimize();
}

/*
 *
 */
mdreal
Array::remove(const mdsize rank) {
  mdreal value = rlnan;
  
  /* Remove value from full version. */
  if(rank < full.size()) {
    value = full[rank];
    this->full[rank] = rlnan;
    if(value != rlnan) this->ndata -= 1;
    if(rank == (full.size() - 1)) {
      full.resize(rank);
      this->nelem -= 1;
    }
  }

  /* Remove value from sparse version. */
  if(sparse.count(rank) > 0) {
    value = sparse[rank];
    (this->sparse).erase(rank);
    this->nelem -= 1;
  }

  /* Optimize data storage. */
  this->optimize();
  return value;
}

/*
 *
 */
mdsize
Array::size() const {
  return ndata;
}

/*
 *
 */
bool
Array::update(const mdsize rank, const mdreal value, const bool flag) {
  if(value == rlnan) return false;
  
  /* Empty array. */
  if((nelem < 1) && (rank < MIN_SIZE_TRIGGER)) {
    full.resize(rank, rlnan);
    full.push_back(value);
    this->ndata = 1;
    this->nelem = full.size();
    return true;
  }

  /* Update non-empty array. */
  if(full.size() > 0) {

    /* Check capacity. */
    if(rank >= nelem) {
      this->nelem = (rank + 1);
      (this->full).resize(nelem, rlnan);
    }

    /* Add a new value. */
    if(full[rank] == rlnan) {
      this->full[rank] = 0.0;
      this->ndata += 1;
    }

    /* Replace or combine with existing value. */
    if(flag) this->full[rank] = value;
    else this->full[rank] += value;
  }
  else {

    /* Check capacity. */
    if(rank >= nelem) this->nelem = (rank + 1);

    /* Check if element exists. */
    if(sparse.count(rank) < 1) {
      sparse[rank] = 0.0;
      this->ndata += 1;
    }
    
    /* Replace or combine with existing value. */
    if(flag) this->sparse[rank] = value;
    else this->sparse[rank] += value;
  }
  return true;
}

/*
 *
 */
vector<mdreal>
Array::values() const {

  /* Return full data as such. */
  if(full.size() > 0) return full;

  /* Expand sparse data to full version. */
  vector<mdreal> values(nelem, rlnan);
  for(map<mdsize, mdreal>::const_iterator it = sparse.begin();
      it != sparse.end(); it++)
    values[it->first] = it->second;
  return values;
}

/*
 *
 */
mdsize
Array::optimize() {

  /* Remove trailing missing elements. */
  while(full.size() > 0) {
    if(full[nelem-1] != rlnan) break;
    this->nelem -= 1;
    full.resize(nelem);
  }

  /* Make sure element count is up-to-date. */
  map<mdsize, mdreal>::reverse_iterator it = sparse.rbegin();
  if(it != sparse.rend()) this->nelem = it->first;

  /* Estimate ratio of data to capacity. */
  mdreal ratio = (ndata + MIN_SIZE_TRIGGER)/(nelem + 1);

  /* Switch to full storage. */
  if((ratio > FULL_RATIO) && (sparse.size() > 0)) {
    this->full = this->values();
    (this->sparse).clear();
  }

  /* Switch to sparse storage. */
  if((ratio < SPARSE_RATIO) && (full.size() > 0)) {
    for(mdsize i = 0; i < full.size(); i++)
      if(full[i] != rlnan) this->sparse[i] = full[i];
    (this->full).clear();
  }
  return nelem;
}
